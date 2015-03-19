function solver = Solver(problem)
    
    
    switch(lower(problem.solver.name))
        
        case('jfnk')
            
            solver.iterationMaximum    = 500    ;
            solver.tolerance           = 1E-10  ;
            solver.gmres.tolerance     = 1E-12  ;
            solver.gmres.epsilon       = 1E-8   ;
            solver.gmres.restart       = -1     ;
            solver.gmres.nu            = 0.15   ;
            solver.backtracker.relax   = 0.5    ;
            
    end
    
    
    if isfield(problem.solver,'preconditioner')
        solver.preconditioner = IntrepidTwilight.executive.Preconditioner(problem);
    else
        solver.preconditioner.apply  = @(x)    x    ;
        solver.preconditioner.update = @(dum) [ ]   ;
    end
            
    if not(isfield(problem.solver,'guard'))
        solver.guard.value = @(x)     x     ;
        solver.guard.step  = @(x,dx) dx     ;
    end
    
end