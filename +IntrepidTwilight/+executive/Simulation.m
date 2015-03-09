function sim = Simulation(problem)
    
    %   Build semidiscretization
    problem.semidiscretization.closure =...
        IntrepidTwilight.executive.build('semidiscretization',problem);
    sim.f = problem.semidiscretization.closure ;
    
    %   Build time stepper
    problem.timeStepper.closure = IntrepidTwilight.executive.build('timestepper',problem);
    sim.ts = problem.timeStepper.closure;
    
    
    %   Build residual
    sim.r = IntrepidTwilight.executive.build('residual',problem);
    
    
    %   Build preconditioner
    sim.pc = IntrepidTwilight.executive.build('preconditioner',problem);
    
    
    %   Add the finalized problem and other fields
    sim.problem = problem;
    sim.run     = @(timeSpan,saveRate,dt) run(timeSpan,saveRate,dt);
    sim.getData = @() getRunData();
    
    %   Initialize closure variables
    tSave = 0;
    qSave = 0;
    
    
    function [] = run(timeSpan,dt,saveRate)
        
        tStart  = timeSpan(1);
        tFinish = timeSpan(2);
        t       = tStart     ;
        
        %   Times of saved data
        tSave = (tStart:saveRate:tFinish)';
        if (tSave ~= tFinish)
            tSave = [tSave ; tFinish];
        end
        
        
        %   Initialize data
        qSave(problem.miscellaneous.nEq,numel(tSave)) = 0   ;
        qSave(:,1) = problem.initialState.q0                ;
        q          = 1.0001*qSave(:,1)                      ;
        
        %   Initialize preconditioner
        pc = sim.pc(dt);
        pc.update(qSave(:,1));
        
        %   Save index
        k = 2;
        
        while (t <= tFinish)
            
            if     ( abs((t+dt) - tSave(k)) > eps() )
                step     = dt       ;
                saveData = false()  ;
            elseif ( abs((t+dt) - tSave(k)) < eps() )
                step     = dt   ;
                saveData = true ;
            else
                step     = tSave(k) - t ;
                saveData = true         ;
            end

            
            %   Calculate next time value
            pc = sim.pc(step);
            pc.update(q);
            [q,stats] = IntrepidTwilight.TenaciousReduction.JFNK(q,sim.r(step),...
                1E-8,problem.constraint,pc);
            t = t + step;
                
            if saveData
                qSave(:,k) = q;
                k          = k + 1;
            end

            sim.ts.qUpdate(q);
            
            fprintf('%5.2E seconds: %3G iterations, %5.2E residual\n',t,stats.iterations,stats.norm);
            
        end
        
    end
    
    function [t,q] = getRunData()
        t = tSave;
        q = qSave;
    end
    
    
end