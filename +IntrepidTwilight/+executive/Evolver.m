function evolver = Evolver(config)
    
    %   Inherit
    evolver = IntrepidTwilight.executive.Component();
    evolver = evolver.changeID(evolver,'evolver','evolver');
    
    %   Set default values
    evolver.set('time.span'          , [0,1]            );
    evolver.set('time.terminator'    , @(varargin) false);
    evolver.set('time.step.maximum'  , 1                );
    evolver.set('time.step.minimum'  , 0                );
    evolver.set('time.step.goal'     , 0.1              );
    evolver.set('initialCondition'   , 0                );
    evolver.set('saveRate'           , 0.1              );
    evolver.set('tolerance.relative' , 1E-3             );
    evolver.set('tolerance.absolute' , 1E-6             );
    %
    %   Overwrite defaults at construction
    if (nargin >= 1)
        evolver.set(config);
    end
    
    
    
    %   Dependencies and Binder
    evolver.set('dependencies',{'state'})   ;
    evolver.bind = @(object) bind(object)   ;
    state        = []                       ;
    %
    function [] = bind(object)
        if isstruct(object)
            switch(lower(object(1).type))
                case('state')
                    state = object;
            end
        end
    end
    
    
    
    
    
    %   Initialize private variables for the times and values
    times = [];
    qs    = [];
    Dqs   = [];
    
    
    
    % ------------------------------------------- %
    %                Public methods               %
    % ------------------------------------------- %
    
    evolver.run     = @() evolve()  ;
    evolver.evolve  = @() evolve()  ;
    function [] = evolve()
        
        %   Pull runtime information
        timeSpan          = evolver.get('time.span')            ;
        terminator        = evolver.get('time.terminator')      ;
        stepSave          = evolver.get('saveRate')             ;
        stepMin           = evolver.get('time.step.minimum')    ;
        stepMax           = evolver.get('time.step.maximum')    ;
        relativeTolerance = evolver.get('tolerance.relative')   ;
        absoluteTolerance = evolver.get('tolerance.absolute')   ;
        fallbacks         = zeros(1,25)                         ;
        dtFall            = zeros(1,8)                          ;
        
        
        %   Deal with dynamic maximum time step
        if not(any(size(stepMax) == [0,0]))
            if isscalar(stepMax)
                stepMaxFun = @(t) stepMax;
            else
                if size(stepMax,1) == 1 %   One entry table
                    stepMaxFun = @(t) stepMax(1,2);
                elseif size(stepMax,1) >= 2
                    stepMaxFun = @(t) interp1(stepMax(:,1),stepMax(:,2),t,'previous',stepMax(end,2));
                end
            end
        else
            error('IntrepidTwilight:executive:Evolver:emptyTimeStepTable',...
                'Time Step Tables must have no empty dimensions');
        end
        


        %   Create vector of save times
        times = (timeSpan(1):stepSave:timeSpan(2)).';
        if (times(end) ~= timeSpan(2))
            times = [times ; timeSpan(2)];
        end
        nSave = numel(times);


        %   Initial data and allocation
        IC  = evolver.get('initialCondition')   ;
        qs  = IC(:,ones(1,nSave))               ;
        Dqs = qs                                ;


        %   Prepare the state for evoltuion
        state.prepare(qs(:,1),timeSpan(1));
        
        
        % ------------------------------------------- %
        %                  Time March                 %
        % ------------------------------------------- %
        %   Take a single step at stepMin to initialize things and
        %   let the dt coast up through the error estimates returned
        %   by State
        t         = times(1)                ;
        q         = qs(:,1)                 ;
        dt        = stepMin                 ;
        [~,Dq]    = state.update(q,t,0)     ;
        iT        = (1:nSave).'             ;
        saved     = [true;false(nSave-1,1)] ;
        stopMarch = terminator(t,q,Dq)      ;
        
        
        while t<=times(end) && not(stopMarch)
            
            while (t <= times(find(not(saved),1,'first'))) && not(stopMarch)
                
                %   Update from last iteration (wasted assignment for first iteration only)
                qold  = q   ;
                Dqold = Dq  ;
                
                %   Solve, adjudicate, and check for termination
                [q,Dq,stats] = state.update(qold,t,dt)                          ;
                [q,Dq,t,dt]  = adjudicateSolution(q,Dq,qold,Dqold,t,dt,stats)   ;
                
                if (norm(q-qold,1) > eps())
                    stopMarch = terminator(t,q,Dq);
                else
                    stopMarch = false;
                end
                
            end
            
            
            %   Determine intermediate values from a Hermite interpolation
            inPast = times <= t;
            iSave  = iT( inPast & not(saved) ).'; % Transposed for for-loop indexing
            if any(iSave)
                qs(:,iT(iSave)) = ...
                    IntrepidTwilight.ConvenientMeans.hermiteInterpolation(...
                    [t-dt;t],[qold,q],[Dqold,Dq],times(iSave));
                
                %   Calculate intermediate derivatives
                for k = iSave
                    [~,Dqs(:,k),~] = state.update(qs(:,k),times(k),0);
                end
                
                %   Contract save indices
                saved = saved | inPast;
            end
            
            
            %   Handle early termination masking
            if stopMarch 
                mask  = times < t;
                times = [times(mask) ; t];
                qs    = [qs(:,mask)  , q];
                Dqs   = [Dqs(:,mask) , Dq];
            end
            
        end




        % ------------------------------------------------------------ %
        %                      Time Step Adjuster                      %
        % ------------------------------------------------------------ %
        function [q,Dq,t,dt] = adjudicateSolution(q,Dq,qold,Dqold,t,dt,stats)
            
            if  (stats.relativeErrorEstimate > relativeTolerance) || ...
                (stats.absoluteErrorEstimate > absoluteTolerance) || ...
                any(strcmpi(stats.returnStatus,...
                    {'TooSmallStepNorm','IterationMaximumReached','HookExitRequest'}))

                % -------------------------------------------------- %
                %                        Failed                      %
                % -------------------------------------------------- %
                %
                %   Adjust time step
                fallbacks = [1,fallbacks(1:24)]  ;
                dt        = max([dt*2^(-sum(fallbacks(1:5))),stepMin]) ;
                dtFall    = dt;
                
                %   Return old values
                q          = qold;
                Dq         = Dqold;
                fprintf('Fallback\n');

            else

                % -------------------------------------------------- %
                %                       Suceeded                     %
                % -------------------------------------------------- %
                %
                %   Update
                t  = t + dt;
                fprintf(...
                    'Converged:: time: %9.7e;  dt: %9.7e;  norm: %9.7e;  iterations: %g\n',...
                    t,dt,stats.norm(end),stats.iterations);
                
                %   Adjust time step
                fallRate  = sum(fallbacks)/numel(fallbacks);
                dt        = min([...
                            2^(1-fallRate) * dt * (fallRate<=0.1) + ...
                            0.5 * dtFall(1)     * (fallRate> 0.1) ,...
                            stepMaxFun(t)]);
                fallbacks = [0,fallbacks(1:24)]  ;
            end
        end

    end
    
    % -------------------------
    
    evolver.getData = @() getData();
    function varargout = getData()
        switch (nargout)
            case 2
                varargout = {qs,times};
            case 3
                varargout = {qs,Dqs,times};
        end
    end
    
    
end