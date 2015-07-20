function jfnk = JFNK(residual,preconditioner)
    
    %   Default parameters
    jfnk.tolerance.residual       = 1E-7   ;
    jfnk.tolerance.stepSize       = 1E-7   ;
    jfnk.maximumIterations        = 100    ;
    jfnk.epsilon                  = 1E-7   ;
    jfnk.gmres.iteration.restarts =  1     ;
    jfnk.gmres.iteration.maximum  = -1     ;
    jfnk.gmres.tolerance          = 1E-10  ;
    jfnk.gmres.nu                 = 0.20   ;
    jfnk.newton.relax.over        = 1.1    ;
    
    %   Hook initializations
    jfnk.hook.presolve  = @(x) [];
    jfnk.hook.postsolve = @(x) [];
    jfnk.hook.prestep   = @(x) [];
    jfnk.hook.poststep  = @(x) [];



    %   Bind on construction if passed
    if (nargin >= 1) && isstruct(residual)
        bind(residual);
    end
    if (nargin >= 2) && isstruct(preconditioner)
        bind(preconditioner);
    end
    
    
    
    %   Closure variable declarations for GMRES
    Z     = 0   ;   % Update's basis vectors
    H     = 0   ;   % Householder vectors for projections
    R     = 0   ;   % Upper-triangular matrix for least squares problem
    Zeros = 0   ;   % Persistent Zero matrix
    alpha = 0   ;   % Vector of projected residuals
    
    
    
    %   Public methods
    jfnk.type  = 'solver'                   ;
    jfnk.is    = @(s) strcmpi(s,jfnk.type)  ;
    jfnk.set   = @(key,value) set(key,value);
    jfnk.bind  = @(object) bind(object)     ;
    jfnk.solve = @(x) solve(x)              ;
    
    
    
    % ================================================================= %
    %                               Allocater                           %
    % ================================================================= %
    arraysAreNotAllocated = true();
    function [] = allocateWorkArrays(x)
        
        %   Determine row count
        nRows = numel(x);
        
        %   Determine column count
        if (jfnk.gmres.iteration.maximum == -1)
            nCols                        = nRows;
            jfnk.gmres.iteration.maximum = nCols;
        else
            nCols = jfnk.gmres.iteration.maximum;
        end
        
        %   Allocate
        Z(nRows,nCols)   = 0        ;
        H(nRows,nCols)   = 0        ;
        R(nRows,nCols)   = 0        ;
        Zeros(nRows-1,1) = 0        ;
        alpha            = [0;Zeros];
        
        %   Flip the switch
        arraysAreNotAllocated = false();
    end
    
    
    
    
    
    % ================================================================= %
    %                       Mutators/Re-binders                         %
    % ================================================================= %
    function [] = bind(object)
        if isstruct(object)
            switch(object(1).type)
                case('residual')
                    residual = object;

                case('preconditioner')
                    preconditioner = object;
            end
        end
    end
    function [] = set(key,value)
        keys = strsplit(key,'.')                ;
        jfnk = setfield(jfnk,{1},keys{:},value) ;
    end
    
    
    
    
    
    % ================================================================= %
    %                    Fully-Coupled Solver                           %
    % ================================================================= %
    
    function [xNL,stats] = solve(xNL)
        
        %   Allocate arrays if not done already
        if arraysAreNotAllocated
            allocateWorkArrays(xNL);
        end
        
        
        % Hook
        
        
        % ----------------------------------------------- %
        %                Non-Linear Iterations            %
        % ----------------------------------------------- %
        
        %   Hook
        jfnk.hook.presolve(xNL);


        %   Allocate stats struct
        stats.iterations                     = 0;
        stats.norm(jfnk.maximumIterations,1) = 0;

        %   Initialize
        [xNL,rNL]     = residual.guard.state(xNL)           ;
        rNLnorm       = norm(rNL,2)                         ;
        stats.norm(1) = rNLnorm                             ;
        notConverged  = rNLnorm > jfnk.tolerance.residual   ;
        preconditioner.initialize(xNL)                      ;
        
        
        %   Iterate
        while notConverged
            
            %   Hook
            jfnk.hook.prestep(xNL);


            %   Evaluate initial residual
            rNL     = residual.value(xNL)   ;
            rNLnorm = norm(rNL,2)           ;
            
            %   Take a step
            [xNL,rNLnorm,dxNorm] = nonlinearStep(xNL,rNL,rNLnorm);
            
            %   Update to new state
            stats.iterations             = stats.iterations + 1 ;
            stats.norm(stats.iterations) = rNLnorm              ;
            preconditioner.update(xNL)                          ;


            %   Display run-time stuff
            Show(rNLnorm);


            %   Iteration clean-up
            normNotDone            =     rNLnorm      >= jfnk.tolerance.residual            ;
            stepNotDone            =      dxNorm      >= jfnk.tolerance.stepSize            ;
            belowMaximumIterations = stats.iterations <= jfnk.maximumIterations             ;
            notConverged           = (normNotDone || stepNotDone) && belowMaximumIterations ;


            % Hook
            jfnk.hook.poststep(xNL);
            
        end
        
        % Hook
        jfnk.hook.postsolve(xNL);
        
        
        %   Contract vector to the number of actuall iterations
        stats.norm = stats.norm(1:stats.iterations);
        
    end





    % ================================================================= %
    %                 Sub-parts of Nonlinear Solvers                    %
    % ================================================================= %
    
   
    function [xNL,rNLnorm,dxNorm] = nonlinearStep(xNLm1,rNLm1,rNLnorm)
        
        %   Advance in a descent direction
        [xNL,rNL,rNLnorm,dxNorm] = quasiNewtonUpdate(xNLm1,rNLm1,rNLnorm);
        
        
        %   Extrapolate solution
        dxExtrap    = rNL .* (xNL - xNLm1)./(rNL - rNLm1);
        xExtrap     = xNL - dxExtrap            ;
        rExtrap     = residual.value(xExtrap)   ;
        rExtrapNorm = norm(rExtrap,2)           ;
        if rExtrapNorm < 0.9*rNLnorm
            xNL     = xExtrap       ;
            rNLnorm = rExtrapNorm   ;
        end

    end
    
    
    
    function [xNew,rNew,rNewNorm,dxNorm] = quasiNewtonUpdate(xOld,rOld,rOldNorm)
        
        % Solve linear system to within linearTolerance
        dx = GMRES(xOld,rOld,rOldNorm);
        
        %   Relax the step size to the gaurded value
        [dx,rNew] = residual.guard.step(xOld,dx);
        
        %   Set-up
        rNewNorm   = norm(rNew,2)               ;
        notReduced = rNewNorm > rOldNorm        ;
        
        
        if notReduced
            
            %   If full step didn't reduce rssidual or produced NaNs,
            %   search for a better, smaller step.
            alphaMin = ...
                IntrepidTwilight.ConvenientMeans.goldenSectionMinimizeLazily(...
                    @(alpha) norm(residual.value(xOld - alpha*dx),2),...
                        [0,1],[rOldNorm,rNewNorm],0.01);
            rNew     = residual.value(xOld - alphaMin*dx)   ;
            rNewNorm = norm(rNew,2)                         ;
            dx       = alphaMin*dx                          ;
            
            
        else
            
            % Over-relaxation
            reduced    = not(notReduced);
            rTrack     = rNew           ;
            rNormTrack = rNewNorm       ;
            while reduced
                dx         = jfnk.newton.relax.over * dx        ;
                rNew       = residual.value(xOld - dx)   ;
                rNewNorm   = norm(rNew,2)                       ;
                reduced    = rNewNorm < 0.9*rNormTrack          ;
                rTrack     = rNew*reduced + rTrack*(1-reduced)  ;
                rNormTrack = norm(rTrack,2)                     ;
            end
            rNew     = rTrack                   ;
            rNewNorm = rNormTrack               ;
            dx       = dx/jfnk.newton.relax.over;
            
        end
        
        % Calculate relaxed x value
        xNew   = xOld - dx      ;
        dxNorm = norm(dx,Inf)   ;
        
    end
    
    
    
    
    
    % ================================================================= %
    %                              GMRES                                %
    % ================================================================= %
    
    %   Outer "Restart" Iterations
    function dx = GMRES(xk,rk,rkNorm)
        dx = 0;
        for k = 1:jfnk.gmres.iteration.restarts
            xk             = xk + dx                ;
            [dx,rk,rkNorm] = GMRESCore(xk,rk,rkNorm);
            
            if rkNorm <= jfnk.gmres.tolerance
                break;
            end
        end
    end
    
    %   Inner Iterations
    function [dx,rk,rkNorm] = GMRESCore(xk,rk0,rk0Norm)
        
        %   Create shortcuts for closure variables
        nu              = jfnk.gmres.nu         ;
        linearTolerance = jfnk.gmres.tolerance  ;
        epsilon         = jfnk.epsilon          ;
        n               = numel(xk)             ;
        I               = 1:n                   ;
        
        
        %   First basis vector for update
        Z(I,1) = rk0 / rk0Norm ;
        
        % First Step (k = 1)
        % Compute J*z1 and store in R
        w      = preconditioner.apply(Z(I,1))                       ;
        R(I,1) = (residual.value(xk + epsilon*w) - rk0) / epsilon   ;
        
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
        
        
        for k = 2:jfnk.gmres.iteration.maximum
            
            % Choose the next basis vector
            if rkNorm <= nu*rkm1Norm
                Z(I,k) = rk/rkNorm  ;   %   GCR (RB-SGMRES) basis
            else
                Z(I,k) = q          ;   %   Simpler GMRES basis
            end
            
            % Compute and store A*zk in R
            w      = preconditioner.apply(Z(I,k))                ;
            R(I,k) = (residual.value(xk + epsilon*w) - rk0) / epsilon   ;
            
            % Apply all previous projections to new the column
            for m = 1:k-1
                R(I,k) = R(I,k) - H(I,m)*(2*H(I,m)'*R(I,k));
            end
            
            % Get the next Householder vector
            h        = R(k:n,k)                     ;
            h        = -Signum(h(1))*norm(h,2)*e - h;
            h        = h ./ norm(h,2)               ;
            H(k:n,k) = h                            ;
            
            %   Apply projection to R to bring it into upper triangular form;
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
        
        if rcond(Rtilde) > 100*eps()
            yk = Rtilde \ alpha(1:k);   % Solve the least-squares problem
        else
            %   Attempt to re-scale
            S  = diag(1./diag(Rtilde))  ;
            yk = (Rtilde*S) \ alpha(1:k);
            yk = S*yk                   ;
        end
        
        dx = Z(I,1:k) * yk                  ;   % Calculate full Newton update
        dx = preconditioner.apply(dx);
    end
end

function s = Signum(s)
    if (s ~= 0)
        s = sign(s);
    else
        s = 1;
    end
end
