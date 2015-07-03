function jfnk = JFNK(residual,preconditioner)
    
    %   Default parameters
    jfnk.tolernace.residual       = 1E-7   ;
    jfnk.tolernace.stepSize       = 1E-7   ;
    jfnk.maximumIterations        = 100    ;
    jfnk.gmres.iteration.restarts = 1      ;
    jfnk.gmres.iteration.maximum  = -1     ;
    jfnk.gmres.tolerance          = 1E-10  ;
    jfnk.gmres.nu                 = 0.20   ;
    jfnk.newton.relax.over        = 1.1    ;
    
    
    
    %   Object bindings
    r  = 0;     rs  = 0;
    pc = 0;     pcs = 0;
    set('residual'      ,residual);
    set('preconditioner',preconditioner);
    
    
    %   Closure variable declarations
    Z     = 0   ;   % Update's basis vectors
    H     = 0   ;   % Householder vectors for projections
    R     = 0   ;   % Upper-triangular matrix for least squares problem
    Zeros = 0   ;   % Persistent Zero matrix
    alpha = 0   ;   % Vector of projected residuals
    
    
    %   Public methods
    jfnk.allocateWorkArrays = @(x) allocateWorkArrays(x)        ;
    jfnk.set                = @(type,object) set(type,object)   ;
    
    if iscell(r)
        jfnk.solve    = @(x) solveAll(x);
        jfnk.solveOne = @(x) solve(x)   ;
    else
        jfnk.solve    = @(x) solve(x);
    end



    % ================================================================= %
    %                               Allocater                           %
    % ================================================================= %
    arraysAreNotAllocated = true();
    function [] = allocateWorkArrays(x)
        
        %   Determine row ...
        if iscell(r)
            nRows = max(cellfun(@(rk,xk) numel(rk(xk)),r,x));
        else
            nRows = numel(r(x));
        end
        %
        %    and then column count
        if (jfnk.gmres.iteration.restart == -1)
            nCols              = nRows;
            jfnk.gmres.maximum = nCols;
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
                if iscell(object);
                    rs = object;
                else
                    r = object;
                end
            case('preconditioner')
                if iscell(object);
                    pcs = object;
                else
                    pc = object;
                end
        end
    end
    
    

    % ================================================================= %
    %                       Non-Linear Solver                           %
    % ================================================================= %
    
    %   Cell intermediary
    function [xNL,stats] = solveAll(xNL)
        n     = length(rs);
        stats = cell(1,n);
        for k = 1:n
            r                 = rs{k};
            pc                = pcs{k};
            [xNL{k},stats{k}] = solve(xNL{k});
        end
    end
    
    %   Core solver
    function [xNL,varargout] = solve(xNL)
        
        %   Parameter defintion/unpacking
        iterMax      = jfnk.iterationMaximum    ;
        rNLTolerance = jfnk.tolerance.residual  ;
        dxTolerance  = jfnk.tolerance.step      ;
        guard        = jfnk.guard               ;
        
        
        %   Allocate arrays if not done already
        if arraysAreNotAllocated
             allocateWorkArrays(xNL);
        end
        
        
        
        
        % ----------------------------------------------- %
        %                    Initialize                   %
        % ----------------------------------------------- %

        %   Let the gaurd make sure xNL is to its liking
        xNL = guard.value(xNL);
        
        %   Initialization
        rNL     = r(xNL)        ;
        rNLnorm = norm(rNL,2)   ;
        pc.initialize(xNL);
        
        % Counter and residual
        rNormVec(iterMax,1) = 0         ;
        rNormVec(1)         = rNLnorm   ;
        iter                = 2         ;
        
        
        
        % ----------------------------------------------- %
        %                Non-Linear Iterations            %
        % ----------------------------------------------- %
        
        %   First iteration
        [dxNorm,xNL,rNL,rNLnorm] = newtonUpdate(xNL,rNL,rNLnorm);
        pc.update(xNL);
        rNormVec(iter) = rNLnorm    ;
        iter           = iter + 1   ;
        normNotDone    = rNLnorm > rNLTolerance;
        stepNotDone    = dxNorm  > dxTolerance ;
        notConverged   = (normNotDone || stepNotDone) && (iter <= iterMax);
        xNLm1          = xNL;
        
        
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
            
            %   Update preconditions
            pc.update(xNL);
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
    end
    
    
    function [dxNorm,xNew,rNew,rNewNorm] = newtonUpdate(xOld,rOld,rOldNorm)
        
        % Solve linear system to within linearTolerance
        dx = GMRES(xOld,rOld,rOldNorm);
        
        %   Relax the step size to the gaurded value
        dx = guard.step(xOld,dx);
        
        
        
        % =============================================== %
        %                       Relaxation                %
        % =============================================== %
        
        %   Set-up
        rNew       = r(xOld - dx)       ;
        rNewNorm   = norm(rNew,2)       ;
        notReduced = rNewNorm > rOldNorm;
        hasNans    = isnan(rNewNorm)    ;
        
        
        
        if notReduced || hasNans
            
            alphaMin = ...
                IntrepidTwilight.ConvenientMeans.goldenSectionSearch(...
                @(alpha) norm(r(xOld - alpha*dx),2),[0,1],[rOldNorm,rNewNorm],0.01);
            rNew     = r(xOld - alphaMin*dx);
            rNewNorm = norm(rNew,2)         ;
            dx       = alphaMin*dx          ;
            
            
        else % Over-relaxation
            
            reduced    = not(notReduced);
            rTrack     = rNew           ;
            rNormTrack = rNewNorm       ;
            while reduced
                dx         = jfnk.newton.relax.over * dx        ;
                rNew       = r(xOld - dx)                       ;
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
        w      = pc.apply(Z(:,1));
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
        alpha(1) = q'*rk0           ;
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
            w      = pc.apply(Z(:,k));
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
        
        if rcond(Rtilde) > 100*eps()
            yk = Rtilde \ alpha(1:k);   % Solve the least-squares problem
        else
            %   Attempt to re-scale
            S  = diag(1./diag(Rtilde))  ;
            yk = (Rtilde*S) \ alpha(1:k);
            yk = S*yk                   ;
        end
        
        dx = Z(:,1:k) * yk  ;   % Calculate full Newton update
        dx = pc.apply(dx)   ;
    end
end

function s = Signum(s)
    if (s ~= 0)
        s = sign(s);
    else
        s = 1;
    end
end
