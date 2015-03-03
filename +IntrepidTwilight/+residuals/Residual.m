function r = Residual(timeStepper)
    
    r = @(q) (q - timeStepper.qLast()) - timeStepper.deltaq(q);
    
end