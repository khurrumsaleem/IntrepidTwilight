function hem = HEM_TwoLevel(hem)
    
    %   Build full object
    spacediscretization = IntrepidTwilight.executive.build(...
        'AdamantWave',hem.discretization.space,hem.get('model'));
    
    
    %   Define CV (outer) block
    CV.sd.update = @(time) spacediscretization.update(time)      ;
    CV.sd.rhs    = @(qCV)  spacediscretization.rhsMassEnergy(qCV);
    %   sdCV.blockDiagonalJacobian
    %   sdCV.jacobian
    
    
    %   Build MC block
    MC.sd.update = @(time) spacediscretization.update(time);
    MC.sd.rhs    = @(qMC)  spacediscretization.rhsMomentum(qMC);
    %   sdMC.blockDiagonalJacobian
    %   sdMC.jacobian
    
    
    
    %   Build time steppers
    CV.ts = IntrepidTwilight.executive.build('TransientStride',hem.discretization.time,CV.sd);
    MC.ts = IntrepidTwilight.executive.build('TransientStride',hem.discretization.time,MC.sd);
    
    
    
    %   Build residuals
    CV.r = IntrepidTwilight.executive.build('executive','Residual',CV.ts);
    MC.r = IntrepidTwilight.executive.build('executive','Residual',MC.ts);
    
    %
    %   Build preconditioners
    CV.pc = IntrepidTwilight.executive.build('executive','Preconditioner',CV.r,'none');
    MC.pc = IntrepidTwilight.executive.build('executive','Preconditioner',MC.r,'none');
    
    
    %   Build solvers
    CV.solver = IntrepidTwilight.executive.build('TenaciousReduction',hem.solver.name,CV.r,CV.pc);
    MC.solver = IntrepidTwilight.executive.build('TenaciousReduction',hem.solver.name,MC.r,MC.pc);
    
    
    %   Build composite residual and solver
    residual.type   = 'residual';
    residual.is     = @(s) strcmpi(s,'residual');
    residual.update = @(v,t,dt) update(v,t,dt);
    
    solver.type  = 'solver';
    solver.is    = @(s) strcmpi(s,'solver');
    solver.solve = @(v) solve(v);
    
    
    %   Build evolver
    evolver = IntrepidTwilight.executive.build('executive','Evolver',solver,residual);
    
    
    %   Set into hem object
    hem.bind(spacediscretization);
    hem.bind(residual)           ;
    hem.bind(solver)             ;
    hem.bind(evolver)            ;
%     hem.bind('timediscretization.CV',CV.ts)              ;
%     hem.bind('timediscretization.MC',MC.ts)              ;
    
    
    
    function [] = update(v,t,dt)
        CV.r.update(v{1},t,dt);
        MC.r.update(v{2},t,dt);
    end
    
    function [value,stats] = solve(v)
        
        %   Current time values
        massEnergy = v{1};
        momentum   = v{2};
        
        %   Set-Up hooks
        CV.solver.set('hook.presolve',@(qCV) presolve());
        CV.solver.set('hook.prestep' ,@(qCV) prestep(qCV));

        %   Solver system
        [massEnergy,stats(1)] = CV.solver.solve(massEnergy);


        %   Final output stuff
        stats.iterations = stats.CV.iterations;
        stats.norm       = stats.CV.notm;
        value = {massEnergy,momentum};


        function [] = presolve()
            spacediscretization.setMass(massEnergy(1:end/2))        ;
            spacediscretization.setEnergy(massEnergy(end/2+1:end))  ;
            spacediscretization.setMomentum(momentum)               ;
            spacediscretization.updateVelocity()                    ;
        end
        
        %   Solve momentum equations
        function [] = prestep(qCV)
            spacediscretization.setMass(qCV(1:end/2))       ;
            spacediscretization.setEnergy(qCV(end/2+1:end)) ;
            spacediscretization.updateThermodynamicState()  ;
            [momentum,stats(2)] = MC.solver.solve(momentum) ;
            spacediscretization.setMomentum(momentum)       ;
            spacediscretization.updateVelocity()            ;
        end
        
    end
    
end