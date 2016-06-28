function gmres = GMRES()
    
    %   Default parameters
    gmres = IntrepidTwilight.TenaciousReduction.Solver();
    gmres = gmres.changeID(gmres,'GMRES');
    
    %   Parameters
    gmres.set('maximumIterations'  ,   100   )  ;
    gmres.set('iterations.restarts',    1    )  ;
    gmres.set('iterations.maximum' ,   -1    )  ;
    gmres.set('tolerance'          , 1.0E-10 )  ;
    gmres.set('nu'                 ,   0.2   )  ;
    gmres.set('nu'                 ,   0.2   )  ;
    gmres.set('MatrixVectorProduct', @(v) v  )  ;
    gmres.set('Preconditioner'     , @(v) v  )  ;
    
    
    
    %   Prepare
    params         = []             ;
    mvProduct      = []             ;
    preconditioner = []             ;
    gmres.prepare  = @() prepare()  ;
    function [] = prepare()
        params         = gmres.get()                ;
        mvProduct      = params.MatrixVectorProduct ;
        preconditioner = params.Preconditioner      ;
    end





    % ================================================================= %
    %                            GMRES (Outer)                          %
    % ================================================================= %
    gmres.solve = @(x0) solve(x0);
    function dx = solve(xk,rk,rkNorm)
        dx = 0;
        for k = 1:params.iterations.restarts
            xk             = xk + dx                ;
            [dx,rk,rkNorm] = GMRESInner(xk,rk,rkNorm);
            
            if rkNorm <= params.tolerance
                break;
            end
        end
    end
    
    
    % ================================================================= %
    %                            GMRES (Inner)                          %
    % ================================================================= %
    function [dx,rk,rkNorm] = GMRESInner(xk,rk0,rk0Norm)
        
        %   Create shortcuts for closure variables
        n               = numel(xk)         ;
        I               = 1:n               ;
        
        
        %   First basis vector for update
        Z(I,1) = rk0 / rk0Norm ;
        
        % First Step (k = 1)
        % Compute J*z1 and store in R
        w      = preconditioner(Z(I,1)) ;
        R(I,1) = mvProduct(w)           ;
        
        % Compute Householder vector to bring R(:,1) into upper triangular form
        e      = [1 ; Zeros(1:n-1)]                 ;
        h      = R(I,1)                             ;
        h      = -Signum(h(1)) * norm(h,2) * e - h  ;
        H(I,1) = h / norm(h,2)                      ;
        
        % Apply projection to R to bring it into upper triangular form
        R(I,1) = R(I,1) - 2 * H(I,1) * (H(I,1)'*R(I,1));
        
        % Get the first column of the unitary matrix
        q = e - 2 * H(I,1) * (H(I,1)'*e);
        e = e(1:n-1);
        
        % Residual update
        alpha(1) = q'*rk0           ;
        rk       = rk0 - alpha(1)*q ;
        
        % Assign residual norms to determine which basis to use
        rkm1Norm = rk0Norm      ;
        rkNorm   = norm(rk,2)   ;
        
        
        for k = 2:params.gmres.iteration.maximum
            
            % Choose the next basis vector
            if rkNorm <= params.nu*rkm1Norm
                Z(I,k) = rk/rkNorm  ;   %   GCR (RB-SGMRES) basis
            else
                Z(I,k) = q          ;   %   Simpler GMRES basis
            end
            
            % Compute and store A*zk in R
            w      = preconditioner(Z(I,k)) ;
            R(I,k) = mvProduct(w)           ;
            
            % Apply all previous projections to new the column
            for m = 1:k-1
                R(I,k) = R(I,k) - H(I,m)*(2*H(I,m)'*R(I,k));
            end
            
            
            % Get the next Householder vector
            h        = R(k:n,k)                     ;
            h        = -Signum(h(1))*norm(h,2)*e - h;
            h        = h ./ norm(h,2)               ;
            H(k:n,k) = h                            ;
            
            
            %   Apply new projection to R to bring it into upper triangular form;
            for m = 1:k
                R(I,m) = R(I,m) - 2 * H(I,k) * (H(I,k)'*R(I,m));
            end
            
            
            % Get the k-th column of the current unitary matrix
            q = [Zeros(1:k-1) ; e - 2*h*(h'*e)];
            for m = k-1:-1:1
                q = q - 2*H(I,m)*(H(I,m)'*q);
            end
            e = e(1:end-1);
            
            % Update residual
            alpha(k) = q'*rk            ;
            rk       = rk - alpha(k)*q  ;
            
            % Update residual norms
            rkm1Norm = rkNorm;
            rkNorm   = norm(rk,2);
            
            % Solve least-squares problem
            if rkNorm < params.tolerance
                break;
            end
            
        end
        
        % Update to x
        Rtilde = triu(R(1:k,1:k));
        
        if rcond(Rtilde) > 100*eps()
            yk = Rtilde \ alpha(1:k);   % Solve the least-squares problem
        else
            %   Attempt to re-scale
            S  = diag(1./diag(Rtilde))  ;
            yk = (Rtilde*S) \ alpha(1:k);
            yk = S*yk                   ;
        end
        
        dx = Z(I,1:k) * yk                  ;   % Calculate full Newton update
        dx = preconditioner(dx);
    end
    
    
    
end


