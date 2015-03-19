function [xNL,varargout] = JFNK(x0,r,epsilon,constraint,preconditioner)

    
    % ================================================================= %
    %                               Set-Up                              %
    % ================================================================= %
    
    % Length
    N    = length(x0)   ;
    Nmax = N            ;
    
    % Constraint check
    if (nargin < 4) || not(isa(constraint,'function_handle'))
        notConstrained = @(x) false();
    else
        notConstrained = @(x) any(not(constraint(x)));
    end
    
    % Tolerances
    LinearTolerance    = 1E-13;
    NonlinearTolerance = 1E-6 ;

    % Matrix allocation
    Z(N,Nmax) = 0   ;   % Update's basis vectors
    H(N,Nmax) = 0   ;   % Householder vectors for projections
    R(N,Nmax) = 0   ;   % Upper-triangular matrix for least squares problem
    
    % Unit vector used for projections
    Zeros(N-1,1) = 0        ;
    alpha        = [0;Zeros];   % Vector of projected residuals

    % Threshold parameter used to determine which basis to use in the GMRES iterations
    nu = 0.15;




    % ================================================================= %
    %                            JFNK Iteration                         %
    % ================================================================= %
    
    % Let x = x0
    xNL     = x0;
    
    %   Relax the step size for a physical solution
    relaxor = 0.999;
    while notConstrained(xNL)
        xNL = relaxor * xNL;
    end


    % Initial r0
    rNL     = -r(xNL)     ;
    rNormNL = norm(rNL,2)  ;
    
    % Determine if the loop is needed
    NotDone          = rNormNL > NonlinearTolerance ;
    stagnantResidual = false()                      ;

    relaxor = 0.5;
    
    
    %   Initialize preconditioner
    preconditioner.update(xNL);
    
    % Counter and residual
    n                     = 1;
    residualVector(100,1) = 0;
    
    while NotDone
        
        % Solve linear system
        dx = GMRES(xNL,rNL,rNormNL);   % Solve the linear system to LinearTolerance


        %   Relax the step size for a physical solution
        while notConstrained(xNL + dx) && (max(abs(dx)) > 1E-12)
            dx = relaxor * dx;
        end


        % Backtracker
        rNLnew     = r(xNL + dx)            ;
        rNormNLnew = norm(rNLnew,2)         ;
        notDone    = rNormNLnew > rNormNL   ;
        while notDone
            dx         = relaxor * dx;
            rNLnew     = r(xNL + dx);
            rNormNLnew = norm(rNLnew,2);
            notDone    = (rNormNLnew > rNormNL) && max(abs(dx)) > 1E-12;
        end
        xNL  = xNL + dx ;   % Calculate relaxed x value

        
        % Check non-linear residual
        rNL               = -rNLnew     ;
        rNormNL           = rNormNLnew  ;
        residualVector(n) = rNormNL     ;
        
        %   Allow preconditioner to do some post-update work
        preconditioner.update(xNL);
        
        
        NotDone = (rNormNL > NonlinearTolerance) && not(stagnantResidual) ;
        n       = n + 1             ;
        fprintf('\t\t\t%5.2E\t%5.2E\n',rNormNL,norm(dx,2));
    end
    
    if (nargout > 1)
        stats.iterations = n - 1 ;
        stats.norm       = rNormNL      ;
        varargout{1}     = stats        ;
    end
    
    
    % ================================================================= %
    %                          GMRES SubFunction                        %
    % ================================================================= %
    function dx = GMRES(xk,rk0,rk0Norm)
        
        Z(:,1) = rk0 / rk0Norm ; % First basis vector for update
        
        % First Step (k = 1)
        % Compute J*z1 and store in R
        w      = preconditioner.apply(Z(:,1));
        R(:,1) = (r(xk + epsilon*w) + rk0) / epsilon;
        
        % Compute Householder vector to bring R(:,1) into upper triangular form
        e      = [1 ; Zeros]                        ;
        h      = R(:,1)                             ;
        h      = -Signum(h(1)) * norm(h,2) * e - h  ;
        H(:,1) = h / norm(h,2)                      ;
        
        % Apply projection to R to bring it into upper triangular form
        R(:,1) = R(:,1) - 2 * H(:,1) * (H(:,1)'*R(:,1));
        
        % Get the first column of the unitary matrix
        q = e - 2 * H(:,1) * (H(:,1)'*e);
        e = e(1:N-1);
        
        % Residual update
        alpha(1) = q'*rk0           ;
        rk       = rk0 - alpha(1)*q ;
        
        % Assign residual norms to determine which basis to use
        rkm1Norm = rk0Norm      ;
        rkNorm   = norm(rk,2)   ;
        
        
        for k = 2:Nmax
            
            % Choose the next basis vector
            if rkNorm <= nu*rkm1Norm
                Z(:,k) = rk/rkNorm  ;   %   GCR (RB-SGMRES) basis
            else
                Z(:,k) = q          ;   %   Simpler GMRES basis
            end
            
            % Compute and store A*zk in R
            w      = preconditioner.apply(Z(:,k));
            R(:,k) = (r(xk + epsilon*w) + rk0) / epsilon;
            
            % Apply all previous projections to new the column
            for m = 1:k-1
                R(:,k) = R(:,k) - H(:,m)*(2*H(:,m)'*R(:,k));
            end
            
            % Get the next Householder vector
            h        = R(k:N,k)                     ;
            h        = -Signum(h(1))*norm(h,2)*e - h;
            h        = h ./ norm(h,2)               ;
            H(k:N,k) = h                            ;
            
            %   Apply projection to R to bring it into upper triangular form;
            for m = 1:k
                R(:,m) = R(:,m) - 2 * H(:,k) * (H(:,k)'*R(:,m));
            end
            
            % Get the k-th column of the current unitary matrix
            q = [Zeros(1:k-1) ; e - 2*h*(h'*e)];
            for m = k-1:-1:1
                q = q - 2*H(:,m)*(H(:,m)'*q);
            end
            e = e(1:end-1);
            
            % Update residual
            alpha(k) = q'*rk            ;
            rk       = rk - alpha(k)*q  ;
            
            % Update residual norms
            rkm1Norm = rkNorm;
            rkNorm   = norm(rk,2);
            
            % Solve least-squares problem
            if rkNorm/norm(xNL,2) < LinearTolerance
                break;
            end
            
        end
        
        % Update to x
        yk = triu(R(1:k,1:k)) \ alpha(1:k)  ;   % Solve the least-squares problem
        dx = Z(:,1:k) * yk                  ;   % Calculate full Newton update
        dx = preconditioner.apply(dx)       ;
    end
end

function s = Signum(s)
    if (s ~= 0)
        s = sign(s);
    else
        s = 1;
    end
end
