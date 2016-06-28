function evolver = Evolver(config)
    
    %   Inherit
    evolver = IntrepidTwilight.executive.Component();
    evolver = evolver.changeID(evolver,'evolver','evolver');
    
    %   Set default values
    evolver.set('time.span'          , [0,1]  );
    evolver.set('time.step.maximum'  , 1      );
    evolver.set('time.step.minimum'  , 0      );
    evolver.set('time.step.goal'     , 0.1    );
    evolver.set('initialCondition'   , 0      );
    evolver.set('saveRate'           , 0.1    );
    evolver.set('tolerance.relative' , 1E-3   );
    evolver.set('tolerance.absolute' , 1E-6   );
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
    saveTimes = []  ;
    qs    = []  ;
    Dqs   = []  ;
    
    
    
    % ------------------------------------------- %
    %                Public methods               %
    % ------------------------------------------- %
    
    evolver.run     = @() evolve()  ;
    evolver.evolve  = @() evolve()  ;
    function [] = evolve()
        
        %   Pull time information
        timeSpan          = evolver.get('time.span')            ;
        stepSave          = evolver.get('saveRate')             ;
        stepMin           = evolver.get('time.step.minimum')    ;
        stepMax           = evolver.get('time.step.maximum')    ;
        relativeTolerance = evolver.get('tolerance.relative')   ;
        absoluteTolerance = evolver.get('tolerance.absolute')   ;
        fallbacks         = zeros(1,25)                         ;
        dtFall            = zeros(1,8)                          ;
        
        %   Create vector of save times
        saveTimes = (timeSpan(1):stepSave:timeSpan(2)).';
        if (saveTimes ~= timeSpan(2))
            saveTimes = [saveTimes ; timeSpan(2)];
        end
        nSave = numel(saveTimes);
        
        %   Initial data and allocation
        IC = evolver.get('initialCondition');
        qs = IC(:,ones(1,nSave));
        q  = qs(:,1);
        
        
        %   Prepare the state for evoltuion
        state.prepare(q,timeSpan(1));
        
        
        % ------------------------------------------- %
        %                  Time March                 %
        % ------------------------------------------- %
        
        %   Take a single step at stepMin to initialize things and
        %   let the dt coast up through the error estimates returned
        %   by State
        t      = timeSpan(1)        ;
        dt     = stepMin            ;
        [~,Dq] = state.update(q,t,0);
        
        
        for k = 2:nSave

            while t <= saveTimes(k)

                %   Update from last iteration (wasted assignment for first iteration only)
                qold  = q   ;
                Dqold = Dq  ;
                
                %   Solve and adjudicate a potential fallback
                [q,Dq,stats] = state.update(qold,t,dt)                          ;
                [q,Dq,t,dt]  = adjudicateSolution(q,Dq,qold,Dqold,t,dt,stats)   ;
                
                if (t > 2.3e-03)
                    g = [];
                end

            end

            qs(:,k) = IntrepidTwilight.ConvenientMeans.hermiteInterpolation(...
                [t-dt;t],[qold,q],[Dqold,Dq],saveTimes(k));
            [~,Dqs(:,k),~] = state.update(qs(:,k),saveTimes(k),0);

        end

        
        function [q,Dq,t,dt] = adjudicateSolution(q,Dq,qold,Dqold,t,dt,stats)
            if  (stats.relativeErrorEstimate > relativeTolerance) || ...
                (stats.absoluteErrorEstimate > absoluteTolerance) || ...
                any(strcmpi(stats.returnStatus,...
                    {'IterationMaximumReached','HookExitRequest'}))
                
                %   Adjust time step
                fallbacks = [1,fallbacks(1:24)]  ;
                dtFall    = dt;
                dt        = max([dt*2^(-sum(fallbacks(1:5))),stepMin]) ;
                
                %   Return old values
                q          = qold;
                Dq         = Dqold;
                fprintf('Fallback\n');

            else
                
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
                            stepMax]);
                fallbacks = [0,fallbacks(1:24)]  ;
            end
        end

    end
    
    % -------------------------
    
    evolver.getData = @() getData();
    function varargout = getData()
        switch (nargout)
            case 2
                varargout = {qs,saveTimes};
            case 3
                varargout = {qs,Dqs,saveTimes};
        end
    end
    
    
end