function simulation = Simulation(solver,residual)

    %   Add the finalized problem and other fields
    simulation.time.span         = [0,1];
    simulation.time.step.maximum = 1    ;
    simulation.time.step.minimum = 0    ;
    simulation.initialCondition  = 0    ;
    simulation.saveRate          = 1    ;
    
    simulation.run     = @(timeSpan,saveRate,dt) run(timeSpan,saveRate,dt);
    simulation.getData = @() getRunData();
    
    %   Initialize closure variables
    times  = 0;
    values = 0;
    
    
    function [] = run(dt)
        
        %   Create time save vector
        tStart   = simulation.time.span(1)  ;
        tFinish  = simulation.time.span(2)  ;
        saveRate = simulation.saveRate      ;
        t        = tStart                   ;
        
        times = (tStart:saveRate:tFinish)';
        if (times ~= tFinish)
            times = [times ; tFinish];
        end
        
        
        %   Initialize data
        values = simulation.initialCondition(:,ones(1,numel(times)));
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
    
    
end