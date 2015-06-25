function solver = Solver(problem)
    
    
    switch(lower(problem.solver.name))
        
        case('jfnk')
            
            solver = problem.solver;
            
            solver.iterationMaximum    = 500    ;
            solver.tolerance.residual  = 1E-6   ;
            solver.tolerance.step      = 1E-8   ;
            solver.gmres.tolerance     = 1E-12  ;
            solver.gmres.epsilon       = 5E-7   ;
            solver.gmres.restart       = -1     ;
            solver.gmres.nu            = 0.15   ;
            solver.newton.relax.under  = 0.5    ;
            solver.newton.relax.over   = 1.5    ;
            
            %   Legacy
            solver.backtracker.relax = solver.newton.relax.under;
            
    end
    
    
    if isfield(problem.solver,'preconditioner')
        solver.preconditioner = IntrepidTwilight.executive.Preconditioner(problem);
    else
        solver.preconditioner = @(dum) struct('apply',@(x)x,'update',@(dum2)[]);
    end
            
    if not(isfield(problem.solver,'guard'))
        solver.guard.value = @(x)     x     ;
        solver.guard.step  = @(x,dx) dx     ;
    end
    
end