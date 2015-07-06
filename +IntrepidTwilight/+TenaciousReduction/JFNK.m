function jfnk = JFNK(residual,preconditioner)
    
    %   Default parameters
    jfnk.tolernace.residual       = 1E-7   ;
    jfnk.tolernace.stepSize       = 1E-7   ;
    jfnk.maximumIterations        = 100    ;
    jfnk.gmres.iteration.restarts =  1     ;
    jfnk.gmres.iteration.maximum  = -1     ;
    jfnk.gmres.tolerance          = 1E-10  ;
    jfnk.gmres.nu                 = 0.20   ;
    jfnk.newton.relax.over        = 1.1    ;
    
    
    %   Containers for object arrays
    residuals       = 0 ;
    preconditioners = 0 ;
    
    
    %   Bind on construction if passed
    if (nargin >= 1) && isstruct(residual)
        set('residual',residual);
        bindSolver();
    end
    if (nargin >= 2) && isstruct(preconditioner)
        set('preconditioner',preconditioner);
    end
    
    
    
    %   Closure variable declarations for GMRES
    Z     = 0   ;   % Update's basis vectors
    H     = 0   ;   % Householder vectors for projections
    R     = 0   ;   % Upper-triangular matrix for least squares problem
    Zeros = 0   ;   % Persistent Zero matrix
    alpha = 0   ;   % Vector of projected residuals
    
    
    
    %   Public methods
    jfnk.allocateWorkArrays = @(x) allocateWorkArrays(x)        ;
    jfnk.set                = @(type,object) set(type,object)   ;
    jfnk.is                 = @(s) strcmpi(s,'solver')          ;
    
    
    
    % ================================================================= %
    %                               Allocater                           %
    % ================================================================= %
    arraysAreNotAllocated = true();
    function [] = allocateWorkArrays(x)
        
        %   Determine row count
        if iscell(x)
            nRows = max(cellfun(@(c) numel(c),x));
        else
            nRows = numel(x);
        end
        
        %   Determine column count
        if (jfnk.gmres.iteration.iteration.maximum == -1)
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
    function [] = set(type,object)
        switch(lower(type))
            case('residual')
                if object.is('residual')
                    if numel(residual) == 1
                        residual = object;
                    else
                        residual  = object(1)   ;
                        residuals = object      ;
                    end
                end
                
            case('preconditioner')
                if object.is('preconditioner')
                    if numel(preconditioner) == 1
                        preconditioner = object;
                    else
                        preconditioner  = object(1) ;
                        preconditioners = object    ;
                    end
                end
        end
    end
    function [] = bindSolver()
        if isstruct(residual) && residual.is('residual')
            if numel(residual) == 1
                jfnk.solve = @(x) solveCoupled(x)   ;
            else
                jfnk.solve = @(x) solveSegregated(x);
            end
        end
    end
    
    
    
    
    
    % ================================================================= %
    %                    Fully-Coupled Solver                           %
    % ================================================================= %
    
    function [xNL,stats] = solveCoupled(xNL)
        
        %   Allocate arrays if not done already
        if arraysAreNotAllocated
            allocateWorkArrays(xNL);
        end
        
        
        % ----------------------------------------------- %
        %                Non-Linear Iterations            %
        % ----------------------------------------------- %
        
        %   Initialize stuff
        [rNL,rNLnorm,stats] = initialize(xNL)                   ;
        notConverged        = rNLnorm > jfnk.tolerance.residual ;
        
        %   Iterate
        while notConverged
            
            %   Take a step
            [xNL,rNL,dxNorm,rNLnorm,stats] = nonlinearStep(xNL,rNL,rNLnorm,stats);
            
            %   Display run-time stuff
            Show(rNLnorm);
            
            %   Iteration clean-up
            normNotDone           =     rNLnorm      >= jfnk.tolerance.residual             ;
            stepNotDone           =      dxNorm      >= jfnk.tolerance.stepSize             ;
            belowIterationMaximum = stats.iterations <= jfnk.iterationMaximum               ;
            notConverged          = (normNotDone || stepNotDone) && belowIterationMaximum   ;
            
        end
        
        
        %   Contract vector to the number of actuall iterations
        stats.rNorm = stats.rNorm(1:stats.iterations)  ;
        
    end
    
    
    
    
    
    % ================================================================= %
    %                       Segregated Solver                           %
    % ================================================================= %
    
    function [xNL,stats] = solveSegregated(xNL)
        
        %   Allocate arrays if not done already
        if arraysAreNotAllocated
            allocateWorkArrays(xNL);
        end
        
        
        % ----------------------------------------------- %
        %                Non-Linear Iterations            %
        % ----------------------------------------------- %
        
        %   Allocate
        n                     = numel(xNL)  ;
        rNL{n,1}              = 0           ;
        rNLnorm(n,1)          = 0           ;
        dxNorm(n,1)           = 0           ;
        stats(n,1).iterations = 0           ;
        
        %   Initialize
        for k = 1:n
            [rNL{k},rNLnorm(k),stats(k)] = initialize(xNL{k});
        end
        notConverged = all(rNLnorm > jfnk.tolerance.residual);
        
        
        %   Iterate
        while notConverged
            
            for k = 1:n
                
                %   Assign this block's objects
                residual       = residuals(k)       ;
                preconditioner = preconditioners(k) ;
                
                
                %   Take a step
                [xNL{k},rNL{k},dxNorm(k),rNLnorm(k),stats(k)] = ...
                    nonlinearStep(xNL{k},rNL{k},rNLnorm(k),stats(k));
                
                %   Display run-time stuff
                fprintf('Block %d:\n\t',k);
                Show(rNLnorm{k});
            end
            
            %   Iteration clean-up
            normNotDone           = all(       rNLnorm       >= jfnk.tolerance.residual)    ;
            stepNotDone           = all(        dxNorm       >= jfnk.tolerance.stepSize)    ;
            belowIterationMaximum =      stats(1).iterations <= jfnk.iterationMaximum       ;
            notConverged          = (normNotDone || stepNotDone) && belowIterationMaximum   ;
            
        end
        
        %   Contract vector to the number of actuall iterations
        for k = 1:n
            stats(k).rNorm = stats(k).rNorm(1:stats(k).iterations);
        end
        
    end
    
    
    
    
    
    % ================================================================= %
    %                 Sub-parts of Nonlinear Solvers                    %
    % ================================================================= %
    
    
    function [rNL,rNLnorm,stats] = initialize(xNL)
        
        %   Guard against the initial value
        [rNL,xNL] = residual.guard.state(xNL);
        
        %   Initialization
        rNLnorm = norm(rNL,2)           ;
        preconditioner.initialize(xNL)  ;
        
        % Counter and residual
        stats.rNorm(jfnk.iterationMaximum,1) = 0        ;
        stats.rNorm(1)                       = rNLnorm  ;
        stats.iterations                     = 1        ;
        
    end
    
    
    
    function [xNL,rNL,rNLnorm,dxNorm,params] = nonlinearStep(xNLm1,rNLm1,rNLnorm,params)
        
        %   Advance in a descent direction
        [xNL,rNL,rNLnorm,dxNorm] = quasiNewtonUpdate(xNLm1,rNLm1,rNLnorm);
        
        
        %   Extrapolate solution
        dxExtrap    = rNL .* (xNL - xNLm1)./(rNL - rNLm1);
        xExtrap     = xNL - dxExtrap            ;
        rExtrap     = residual.value(xExtrap)   ;
        rExtrapNorm = norm(rExtrap,2)           ;
        if rExtrapNorm < 0.9*rNLnorm
            xNL     = xExtrap       ;
            rNL     = rExtrap       ;
            rNLnorm = rExtrapNorm   ;
        end
        
        %   Update preconditions
        preconditioner.update(xNL);
        Show(rNLnorm);
        
        %   Iteration clean-up
        params.iterations  = params.iterations + 1  ;
        params.rNorm(iter) = rNLnorm                ;
    end
    
    
    
    function [xNew,rNew,rNewNorm,dxNorm] = quasiNewtonUpdate(xOld,rOld,rOldNorm)
        
        % Solve linear system to within linearTolerance
        dx = GMRES(xOld,rOld,rOldNorm);
        
        %   Relax the step size to the gaurded value
        [rNew,dx] = residual.guard.step(xOld,dx);
        
        %   Set-up
        rNewNorm   = norm(rNew,2)               ;
        notReduced = rNewNorm > rOldNorm        ;


        if notReduced
            
            %   If full step didn't reduce rssidual or produced NaNs,
            %   search for a better, smaller step.
            alphaMin = ...
                IntrepidTwilight.ConvenientMeans.goldenSectionMinimizeLazily(...
                @(alpha) norm(residual.value(xOld - alpha*dx),2),[0,1],[rOldNorm,rNewNorm],0.01);
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
                rNew       = residual.value(xOld - dx)          ;
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
            
            if rkNorm <= jfnk.gmres.tolernace
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
        
        
        %   First basis vector for update
        Z(:,1) = rk0 / rk0Norm ;
        
        % First Step (k = 1)
        % Compute J*z1 and store in R
        w      = preconditioner.apply(Z(:,1))                       ;
        R(:,1) = (residual.value(xk + epsilon*w) - rk0) / epsilon   ;
        
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
        alpha(1) = q'*rk0           ;
        rk       = rk0 - alpha(1)*q ;
        
        % Assign residual norms to determine which basis to use
        rkm1Norm = rk0Norm      ;
        rkNorm   = norm(rk,2)   ;
        
        
        for k = 2:nMax
            
            % Choose the next basis vector
            if rkNorm <= nu*rkm1Norm
                Z(:,k) = rk/rkNorm  ;   %   GCR (RB-SGMRES) basis
            else
                Z(:,k) = q          ;   %   Simpler GMRES basis
            end
            
            % Compute and store A*zk in R
            w      = preconditioner.apply(Z(:,k))                       ;
            R(:,k) = (residual.value(xk + epsilon*w) - rk0) / epsilon   ;
            
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
        
        if rcond(Rtilde) > 100*eps()
            yk = Rtilde \ alpha(1:k);   % Solve the least-squares problem
        else
            %   Attempt to re-scale
            S  = diag(1./diag(Rtilde))  ;
            yk = (Rtilde*S) \ alpha(1:k);
            yk = S*yk                   ;
        end
        
        dx = Z(:,1:k) * yk              ;   % Calculate full Newton update
        dx = preconditioner.apply(dx)   ;
    end
end

function s = Signum(s)
    if (s ~= 0)
        s = sign(s);
    else
        s = 1;
    end
end
