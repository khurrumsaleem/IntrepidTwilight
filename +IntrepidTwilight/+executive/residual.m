function r = Residual(problem)

    r  = @(dt) @(q) value(q,dt) ;
    
    function r = value(q,dt)
        r = (q - problem.timeStepper.closure.qLast()) - problem.timeStepper.closure.deltaq(q,dt);
    end
    
end