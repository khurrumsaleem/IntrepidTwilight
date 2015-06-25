function simulation = Simulation(problem)
    
    %     sim.problem = problem;
    %
    %
    %     %   Build semidiscretization
    %     sim.f = IntrepidTwilight.executive.build('semidiscretization',problem);
    %     problem.meta__.builtComponents.f = sim.f;
    
    
    
    
    
    % ============================================================= %
    %                       Create Closure                          %
    % ============================================================= %
    
    %   Methods
    simulation.semidiscretization = @(sd) set('semidiscretization',sd);
    simulation.timestepper        = @(ts) set('timestepper',ts);
    simulation.model              = @(md) set('model',md);
    simulation.run                = [];
    simulation.plot               = [];
    
    %   Public properties
    simulation.stats              = [];
    simulation.results            = [];
    
    
    
    %   Private data
    semidiscretization.isSet = false();
    semidiscretization.name  = '';
    semidiscretization.value = [];
    timestepper.isSet        = false();
    timestepper.name         = '';
    timestepper.value        = [];
    model.isSet              = false();
    model.name               = '';
    model.value              = [];
    readyToRun               = semidiscretization.isSet & timestepper.isSet & model.isSet;
    
    
    
    % ============================================================= %
    %                    Initialize Closure                         %
    % ============================================================= %
    
    if (nargin >= 1) && not(isempty(problem))
        for field = fieldnames(problem)
            simulation.(field)(problem.(field));
        end
    end
    
    
    
    
    
    
    function [] = set(component,value)
        
        
        switch(lower(component))
            
            case('semidiscretization')
                % ================================================= %
                %          Semidiscretization Defintion             %
                % ================================================= %
                semidiscretization.name = value;
                
                if model.isSet
                    semidiscretization.value = ...
                        IntrepidTwilight.executive.makeSemidiscretization(semidiscretization.name,model.value);
                    semidiscretization.isSet = true();
                end
                
                
                
                
                
            case('timestepper')
                % ================================================= %
                %                TimeStepper Defition               %
                % ================================================= %
                
                %   Set name
                timestepper.name = value;
                
                
                %   Check if the semidiscretization can and should be built
                if  model.isSet && not(isempty(semidiscretization.name))
                    set('semidiscretization',semidiscretization.name);
                end
                
                
                %   Check if the timestepper can and should be built
                if semidiscretization.isSet
                    
                    if which(['IntrepidTwilight.TransientStride.',timestepper.name])
                        timestepper.value = ...
                            IntrepidTwilight.TransientStride.(timestepper.name)(semidiscretization.value);
                        timestepper.isSet = true;
                    else
                        error('IntrepidTwilight:executive:build:unknownTimestepper',...
                            'The requested timestepper ''%s'' could not be found.',stepper);
                    end
                    
                end
                
                
                
                
                
            case('model')
                % ================================================= %
                %                    Model Defition                 %
                % ================================================= %
                
                %   Define the model
                if isstruct(value) || isa(value,'Model');
                    model.value = value;
                    if isfield(value,'name') || isprop(value,'name')
                        model.name = model.value.name;
                    end
                    model.isSet = true;
                end
                
        end
    end
    
end


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
    q          = qSave(:,1)                             ;
    
    %   Initialize preconditioner
    %         pc = sim.pc(dt);
    %         pc.update(qSave(:,1));
    sim.f.updateTime(tStart);
    
    %   Save index
    k = 2;
    
    while ( k <= numel(tSave) )
        
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
        %             pc = sim.pc(step);
        %             pc.update(q);
        t = t + step;
        sim.f.updateTime(t);
        [q,stats] = IntrepidTwilight.TenaciousReduction.JFNK(1.0001*q,sim.r(step),...
            sim.solver);
        
        if (t > 0.540)
            g = [];
        end
        
        if saveData
            qSave(:,k) = q;
            k          = k + 1;
        end
        
        sim.ts.qUpdate(q);
        fprintf('%5.2E seconds: %3G iterations, %5.2E residual\n',t,stats.iterations,stats.norm(end));
        
    end
    
end



