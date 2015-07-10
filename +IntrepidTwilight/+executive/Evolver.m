function evolver = Evolver(Solver,Residual)
    
    
    %   Bind at construction
    if (nargin >= 1) && not(isempty(Solver))
        set('solver',Solver);
    else
        solver = 0;
    end
    if (nargin >= 2) && not(isempty(Residual))
        set('residual',Residual);
    else
        residual = 0;
    end
    

    %   Add the finalized problem and other fields
    evolver.time.span         = [0,1];
    evolver.time.step.maximum = 1    ;
    evolver.time.step.minimum = 0    ;
    evolver.time.step.goal    = 0.1  ;
    evolver.initialCondition  = 0    ;
    evolver.saveRate          = 1    ;
    
    evolver.set     = @(type,object) set(type,object);
    evolver.run     = @() run();
    evolver.getData = @() getRunData();
    
    %   Initialize closure variables
    times  = 0;
    values = 0;
    
    
    function [] = run()
        
        %   Create time save vector
        tStart   = evolver.time.span(1)     ;
        tFinish  = evolver.time.span(2)     ;
        saveRate = evolver.saveRate         ;
        dt       = evolver.time.step.goal   ;
        t        = tStart                   ;
        
        times = (tStart:saveRate:tFinish)';
        if (times ~= tFinish)
            times = [times ; tFinish];
        end
        
        
        %   Initialize data
        values = evolver.initialCondition(:,ones(1,numel(times)));
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
                k          = k + 1;
            end

            %   Print statistics
            fprintf('%5.2E seconds: %3G iterations, %5.2E residual\n',t,stats.iterations,stats.norm(end));
            
            
        end
        
    end
    
    function [t,q] = getRunData()
        t = times;
        q = values;
    end
    
    
    
    
    %   Late binder
    function [] = set(type,object)
        switch(lower(type))
            case('residual')
                if object.is('residual')
                    residual = object;
                end

            case('solver')
                if object.is('solver')
                    solver = object;
                end
        end
    end
    
    
end