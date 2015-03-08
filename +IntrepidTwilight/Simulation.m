function sim = Simulation(problem)

    sim.problem = problem;


    %   Build semidiscretization
    sim.f = IntrepidTwilight.executive.build('semidiscretization',problem);
    problem.meta__.builtComponents.f = sim.f;
    

end