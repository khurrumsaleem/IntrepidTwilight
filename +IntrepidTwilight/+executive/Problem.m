function p = Problem()

    p.geometry              = [];
    p.initialState          = [];
    p.miscellaneous         = [];
    p.timeStepper.name      = '';
    p.timeStepper.stepSize  = [];
    p.semidiscretization    = '';
    p.solver.name           = '';
    p.solver.preconditioner = '';
    p.solver.guard          = '';
    
end