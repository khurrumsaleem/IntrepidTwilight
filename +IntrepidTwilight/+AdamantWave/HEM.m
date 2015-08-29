function hem = HEM()
    
    %   Inherit
    parent = IntrepidTwilight.executive.Simulation();
    hem    = parent;
    
    
%     %   Set Defaults
%     hem.choose.spacediscretization.Quasi2DUpwind();
%     hem.choose.timediscretization.ImplicitEuler();
%     hem.choose.residual.residual();
%     hem.choose.preconditioner.Preconditioner();
%     hem.choose.solver.JFNK();
%     hem.choose.evolver.evolver();
%     hem.choose.scheme.twolevel();
    
    %   Bind basic AdamantWave module
    hem.bind(IntrepidTwilight.AdamantWave.BasicModel());
    
    
end

function hem = buildTwoLevel(hem)

    %   Build full SD and subassign
    comp = hem.get('spacediscretization');
    sd   = IntrepidTwilight.executive.build(comp.module,comp.name,comp);
        sd.bind(hem.get('model'));
    %
    outer.sd.update = @(time) sd.update(time)      ;
    outer.sd.rhs    = @(qCV)  sd.rhsMassEnergy(qCV);
    %   outer.sd.blockDiagonalJacobian
    %   outer.sd.jacobian
    %
    inner.sd.update = @(time) sd.update(time)       ;
    inner.sd.rhs    = @(qMC)  sd.rhsMomentum(qMC)   ;
    %   inner.sd.blockDiagonalJacobian
    %   inner.sd.jacobian
      
    
    %   Build time steppers
    comp     = hem.get('timediscretization');
    outer.ts = IntrepidTwilight.executive.build(comp.module,comp.name,comp);
        outer.ts.bind(outer.sd);
    inner.ts = IntrepidTwilight.executive.build(comp.module,comp.name,comp);
        inner.ts.bind(inner.sd);
    
    
    
    %   Build residuals
    comp    = hem.get('residual');
    outer.r = IntrepidTwilight.executive.build(comp.module,comp.name,comp);
        outer.r.bind(outer.ts);
    inner.r = IntrepidTwilight.executive.build(comp.module,comp.name,comp);
        inner.r.bind(inner.ts);
    
    %
    %   Build preconditioners
    comp     = hem.get('preconditioner');
    outer.pc = IntrepidTwilight.executive.build(comp.module,comp.name,comp);
        outer.pc.bind(outer.r);
    inner.pc = IntrepidTwilight.executive.build(comp.module,comp.name,comp);
        inner.pc.bind(inner.r);
    
    
    %   Build solvers
    comp         = hem.get('solver');
    outer.solver = IntrepidTwilight.executive.build(comp.module,comp.name,comp);
        outer.solver.bond(outer.r);
        outer.solver.bond(outer.pc);
    inner.solver = IntrepidTwilight.executive.build(comp.module,comp.name,comp);
        inner.solver.bond(inner.r);
        inner.solver.bond(inner.pc);
    
    
    %   Build composite residual and solver
    residual        = IntrepidTwilight.executive.Residual() ;
    residual.update = @(v,t,dt) update(v,t,dt)              ;
    
    solver       = IntrepidTwilight.TenaciousReduction.Solver();
    solver.solve = @(v) solve(v);
    
    
    %   Build evolver
    comp    = hem.get('solver');
    evolver = IntrepidTwilight.executive.build(comp.module,comp.name,comp);
        evolver.bind(solver)    ;
        evolver.bind(residual)  ;
    
    
    %   Set into hem object
    hem.bind(sd);
    hem.bind(residual)           ;
    hem.bind(solver)             ;
    hem.bind(evolver)            ;
%     hem.bind('timediscretization.CV',CV.ts)              ;
%     hem.bind('timediscretization.MC',MC.ts)              ;
    
    
    
    function [] = update(v,t,dt)
        outer.r.update(v{1},t,dt);
        inner.r.update(v{2},t,dt);
    end
    
    function [value,stats] = solve(v)
        
        %   Current time values
        massEnergy = v{1};
        momentum   = v{2};
        
        %   Set-Up hooks
        outer.solver.set('hook.presolve',@(qCV) presolve());
        outer.solver.set('hook.prestep' ,@(qCV) prestep(qCV));

        %   Solver system
        [massEnergy,stats(1)] = outer.solver.solve(massEnergy);


        %   Final output stuff
        stats.iterations = stats.CV.iterations;
        stats.norm       = stats.CV.notm;
        value = {massEnergy,momentum};


        function [] = presolve()
            sd.setMass(massEnergy(1:end/2))        ;
            sd.setEnergy(massEnergy(end/2+1:end))  ;
            sd.setMomentum(momentum)               ;
            sd.updateVelocity()                    ;
        end
        
        %   Solve momentum equations
        function [] = prestep(qCV)
            sd.setMass(qCV(1:end/2))       ;
            sd.setEnergy(qCV(end/2+1:end)) ;
            sd.updateThermodynamicState()  ;
            [momentum,stats(2)] = inner.solver.solve(momentum) ;
            sd.setMomentum(momentum)       ;
            sd.updateVelocity()            ;
        end
        
    end
    
end
