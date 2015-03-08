function item = build(buildWhat,problem)
    
    switch(lower(buildWhat))
        case('semidiscretization')
            item = IntrepidTwilight.executive.makeSemidiscretization(problem);

        case('timestepper')
            
            stepper = problem.timeStepper.name;
            
            if which(['IntrepidTwilight.TransientStride.',stepper])
                f    = problem.semidiscretization.closure;
                dt   = problem.timeStepper.stepSize;
                item = IntrepidTwilight.TransientStride.(stepper)(f,dt);
            else
                error('IntrepidTwilight:executive:build:unknownTimestepper',...
                    'The requested timestepper ''%s'' could not be found.',stepper);
            end
            
            
            
        case('residual')
            
        otherwise
    end
    
end