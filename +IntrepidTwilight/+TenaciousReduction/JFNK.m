function [xNL,varargout] = JFNK(x0,r,parameters)
    
    
    % ================================================================= %
    %                               Set-Up                              %
    % ================================================================= %
    
    %   Parameter defintion/unpacking
    n               = length(x0)                     ;
    iterMax         = parameters.iterationMaximum    ;
    rNLTolerance    = parameters.tolerance.residual  ;
    dxTolerance     = parameters.tolerance.step      ;
    nMax            = parameters.gmres.restart       ;
    nu              = parameters.gmres.nu            ;
    linearTolerance = parameters.gmres.tolerance     ;
    epsilon         = parameters.gmres.epsilon       ;
    gammaUnder      = parameters.newton.relax.under  ;
    gammaOver       = parameters.newton.relax.over   ;
    preconditioner  = parameters.preconditioner      ;
    guard           = parameters.guard               ;
    
    if (nMax == -1)
        nMax = n;
    end
    
    
    % Matrix allocation
    Z(n,nMax) = 0   ;   % Update's basis vectors
    H(n,nMax) = 0   ;   % Householder vectors for projections
    R(n,nMax) = 0   ;   % Upper-triangular matrix for least squares problem
    
    % Unit vector used for projections
    Zeros(n-1,1) = 0        ;
    alpha        = [0;Zeros];   % Vector of projected residuals
    
    % ================================================================= %
    %                            JFNK Iteration                         %
    % ================================================================= %
    
    % Let x = x0
    xNL = x0;
    
    %   Let the gaurd make sure xNL is to its liking
    xNL = guard.value(xNL);
    
    % Initial r0
    rNL     = r(xNL)        ;
    rNLnorm = norm(rNL,2)   ;
    
    %   Initialize preconditioner
    preconditioner.update(xNL);
    
    % Counter and residual
    rNormVec(iterMax,1) = 0         ;
    rNormVec(1)         = rNLnorm   ;
    iter                = 2         ;
    
    
    
    % =========================================================================== %
    %                            Non-linear iterations                            %
    % =========================================================================== %
    
    %   First iteration
    [dxNorm,xNL,rNL,rNLnorm] = newtonUpdate(xNL,rNL,rNLnorm);
    preconditioner.update(xNL);
    rNormVec(iter) = rNLnorm    ;
    iter           = iter + 1   ;
    normNotDone    = rNLnorm > rNLTolerance;
    stepNotDone    = dxNorm  > dxTolerance ;
    notConverged   = (normNotDone || stepNotDone) && (iter <= iterMax);
    xNLm1          = xNL;
    
    
    %     beta = 2;
    
    while notConverged
        %   Get the Newton update with back-tracking
        [dxNorm,xNL,rNL,rNLnorm] = newtonUpdate(xNL,rNL,rNLnorm);
        
        
        %   Extrapolate solution
        dxExtrap    = rNL .* (xNL - xNLm1)./(rNL - rNormVec(iter-1));
        xExtrap     = xNL - dxExtrap    ;
        rExtrap     = r(xExtrap)        ;
        rExtrapNorm = norm(rExtrap,2)   ;
        if rExtrapNorm < 0.9*rNLnorm
            xNL     = xExtrap       ;
            rNL     = rExtrap       ;
            rNLnorm = rExtrapNorm   ;
        end
        
        %   Allow preconditioner to do some post-newton updating
        preconditioner.update(xNL);
        
        Show(rNLnorm);
        
        %   Iteration clean-up
        xNLm1          = xNL        ;
        rNormVec(iter) = rNLnorm    ;
        iter           = iter + 1   ;
        normNotDone    = rNLnorm > rNLTolerance;
        stepNotDone    = dxNorm  > dxTolerance ;
        notConverged   = (normNotDone || stepNotDone) && (iter <= iterMax);
    end
    
    
    if (nargout > 1)
        stats.iterations = iter - 1                     ;
        stats.norm       = rNormVec(1:stats.iterations) ;
        varargout{1}     = stats                        ;
    end
    
    
    
    function [dxNorm,xNew,rNew,rNewNorm] = newtonUpdate(xOld,rOld,rNormOld)
        
        % Solve linear system to within linearTolerance
        dx = GMRES(xOld,rOld,rNormOld);
        
        %   Relax the step size to the gaurded value
        dx = guard.step(xOld,dx);
        
        
        
        % =============================================== %
        %                       Relaxation                %
        % =============================================== %
        
        %   Set-up
        rNew       = r(xOld - dx)         ;
        rNewNorm   = norm(rNew,2)         ;
        notReduced = rNewNorm > rNormOld  ;
        
        if notReduced   % Under-relaxation
            
            while notReduced
                dx         = gammaUnder*dx          ;
                rNew       = r(xOld - dx)           ;
                rNewNorm   = norm(rNew,2)           ;
                notReduced = rNewNorm > rNormOld    ;
            end
            
        else % Over-relaxation

            reduced    = not(notReduced);
            rTrack     = rNew           ;
            rNormTrack = rNewNorm       ;
            while reduced
                dx         = gammaOver * dx                     ;
                rNew       = r(xOld - dx)                       ;
                rNewNorm   = norm(rNew,2)                       ;
                reduced    = rNewNorm < 0.9*rNormTrack          ;
                rTrack     = rNew*reduced + rTrack*(1-reduced)  ;
                rNormTrack = norm(rTrack,2)                     ;
            end
            rNew     = rTrack       ;
            rNewNorm = rNormTrack   ;
            dx       = dx/gammaOver ;

        end

        % Calculate relaxed x value
        xNew   = xOld - dx      ;
        dxNorm = norm(dx,Inf)   ;

    end
    
    
    
    
    
    
    
    
    
    % ================================================================= %
    %                          GMRES SubFunction                        %
    % ================================================================= %
    function dx = GMRES(xk,rk0,rk0Norm)
        
        Z(:,1) = rk0 / rk0Norm ; % First basis vector for update
        
        % First Step (k = 1)
        % Compute J*z1 and store in R
        w      = preconditioner.apply(Z(:,1));
        R(:,1) = (r(xk + epsilon*w) - rk0) / epsilon;
        
        % Compute Householder vector to bring R(:,1) into upper triangular form
        e      = [1 ; Zeros]                        ;
        h      = R(:,1)                             ;
        h      = -Signum(h(1)) * norm(h,2) * e - h  ;
        H(:,1) = h / norm(h,2)                      ;
        
        % Apply projection to R to bring it into upper triangular form
        R(:,1) = R(:,1) - 2 * H(:,1) * (H(:,1)'*R(:,1));
        
        % Get the first column of the unitary matrix
        q = e - 2 * H(:,1) * (H(:,1)'*e);
        e = e(1:n-1);
        
        % Residual update
        alpha(1) = q'*rk0          ;
        rk       = rk0 - alpha(1)*q ;
        
        % Assign residual norms to determine which basis to use
        rkm1Norm = rk0Norm      ;
        rkNorm   = norm(rk,2)   ;
        
        
        for k = 2:nMax
            
            % Choose the next basis vector
            if rkNorm <= nu*rkm1Norm
                Z(:,k) = rk/rkNorm ;   %   GCR (RB-SGMRES) basis
            else
                Z(:,k) = q          ;   %   Simpler GMRES basis
            end
            
            % Compute and store A*zk in R
            w      = preconditioner.apply(Z(:,k));
            R(:,k) = (r(xk + epsilon*w) - rk0) / epsilon;
            
            % Apply all previous projections to new the column
            for m = 1:k-1
                R(:,k) = R(:,k) - H(:,m)*(2*H(:,m)'*R(:,k));
            end
            
            % Get the next Householder vector
            h        = R(k:n,k)                     ;
            h        = -Signum(h(1))*norm(h,2)*e - h;
            h        = h ./ norm(h,2)               ;
            H(k:n,k) = h                            ;
            
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
            alpha(k) = q'*rk           ;
            rk       = rk - alpha(k)*q  ;
            
            % Update residual norms
            rkm1Norm = rkNorm;
            rkNorm   = norm(rk,2);
            
            % Solve least-squares problem
            if rkNorm < linearTolerance
                break;
            end
            
        end
        
        % Update to x
        Rtilde = triu(R(1:k,1:k));
        
        if rcond(Rtilde) > 10*eps()
            yk = Rtilde \ alpha(1:k);   % Solve the least-squares problem
        else
            S  = diag(1./diag(Rtilde));
            yk = (Rtilde*S) \ alpha(1:k);
            yk = S*yk;
        end
        
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
