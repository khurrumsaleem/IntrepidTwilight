function evolver = makeEvolver(sd,alpha)

    %   Instantiate
    ts      = IntrepidTwilight.TransientStride.TRBDF2()     ;
        ts.set('alpha',alpha);
    r       = IntrepidTwilight.executive.Residual()         ;
    pc      = IntrepidTwilight.executive.Preconditioner()   ;
    solver  = IntrepidTwilight.TenaciousReduction.JFNK()    ;
    state   = IntrepidTwilight.executive.State()            ;
    evolver = IntrepidTwilight.executive.Evolver()          ;
    
    %   Bind some stuff
    ts.bind(sd)         ;
    r.bind(ts)          ;
    pc.bind(r)          ;
    solver.bind(r,pc)   ;
    state.bind(r,solver);


    %   Create buffer between the evolver and the state to account for the 
    %   additional Runge-Kutta stage values
    state   = IntrepidTwilight.executive.State();
    state.bind(r,solver);
    inbetween = IntrepidTwilight.executive.State();
    inbetween.prepare = @(q,t) state.prepare(q,t);
    inbetween.update  = @(q,t,dt) update(q,t,dt);
    function [q,t,stats] = update(q,t,dt)
        [q,t,stats] = state.update([q;q],t,dt);
        q = q(1:end/2);
    end
    
    %   Bind it
    evolver.bind(inbetween);    
    
end
