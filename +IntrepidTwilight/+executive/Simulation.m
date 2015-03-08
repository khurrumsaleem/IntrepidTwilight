function sim = Simulation(problem)

    %   Build semidiscretization
    problem.semidiscretization.closure =...
        IntrepidTwilight.executive.build('semidiscretization',problem);
    sim.f = problem.semidiscretization.closure ;
    
    %   Build time stepper
    problem.closures.ts = IntrepidTwilight.executive.build('timestepper',problem);
    sim.ts = problem.closures.ts;


    %   Build residual
    sim.r = IntrepidTwilight.executive.build('residual',problem);
    
    sim.problem = problem;
    
end