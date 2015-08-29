function evolver = Evolver(config)

    %   Inherit
    evolver = IntrepidTwilight.executive.Component();
    evolver = evolver.changeID(evolver,'evolver','evolver');

    %   Set default values
    evolver.set('time.span'        , [0,1]  );
    evolver.set('time.step.maximum', 1      );
    evolver.set('time.step.minimum', 0      );
    evolver.set('time.step.goal'   , 0.1    );
    evolver.set('initialCondition' , 0      );
    evolver.set('saveRate'         , 0.1    );
    %
    %   Overwrite defaults at construction
    if (nargin >= 1)
        evovler.set(config);
    end
    
    

    %   Dependencies and Binder
    evolver.dependencies = {'residual','solver'}    ;
    evolver.bind         = @(object) bind(object)   ;
    %
    function [] = bind(object)
        if isstruct(object)
            switch(object(1).type)
                case('residual')
                    residual = object;
                case('solver')
                    solver = object;
            end
        end
    end


    %   Public methods
    evolver.run     = @() evolve()  ;
    evolver.evolve  = @() evolve()  ;
    evolver.getData = @() getData() ;
    
    %   Initialize private variables
    residual = []   ;
    solver   = []   ;
    times    = []   ;
    values   = []   ;
    
    
    function [] = evolve()
        
        %   Create time save vector
        tSpan    = evolver.get('time.span')     ;
        saveRate = evolver.get('saveRate')      ;
        dt       = evolver.get('time.step.goal');
        t        = tSpan(1)                     ;
        
        times = (tSpan(1):saveRate:tSpan(2))';
        if (times ~= tSpan(2))
            times = [times ; tSpan(2)];
        end
        
        
        %   Initialize data
        IC     = evolver.get('initialCondition');
        values = IC(:,ones(1,numel(times)));
        value  = values(:,1);
        
        %   Save index
        k = 2;
        
        while ( k <= numel(times) )
            
            %    Adjust time step to coincide with write times
            if     ( abs((t+dt) - times(k)) > eps() )
                step     = dt       ;
                saveData = false()  ;
            elseif ( abs((t+dt) - times(k)) < eps() )
                step     = dt   ;
                saveData = true ;
            else
                step     = times(k) - t ;
                saveData = true         ;
            end



            %   Update time
            t = t + step;
            residual.update(value,t,dt);
            
            %   Solve
            [value,stats] = solver.solve(value);



            %   Store
            if saveData
                values(:,k) = q;
                k           = k + 1;
            end

            %   Print statistics
            fprintf('%5.2E seconds: %3G iterations, %5.2E residual\n',t,stats(1).iterations,stats(1).norm(end));
            
            
        end
        
    end
    
    function [t,q] = getData()
        t = times;
        q = values;
    end
    
    
end