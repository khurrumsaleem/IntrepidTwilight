function item = build(buildWhat,problem)
    
    switch(lower(buildWhat))
        case('semidiscretization')
            item = IntrepidTwilight.executive.makeSemidiscretization(problem);

        case('timestepper')
            
            stepper = problem.timeStepper.name;
            
            if which(['IntrepidTwilight.TransientStride.',stepper])
                item = IntrepidTwilight.TransientStride.(stepper)(problem);
            else
                error('IntrepidTwilight:executive:build:unknownTimestepper',...
                    'The requested timestepper ''%s'' could not be found.',stepper);
            end
            
            
            
        case('residual')
            item = IntrepidTwilight.executive.Residual(problem);
            
            
        case('solver')
            item = IntrepidTwilight.executive.Solver(problem);
            
            
            
        otherwise
    end
    
end