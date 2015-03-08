function r = residual(problem)

    %   Build semi-discretization closure
    closure = ['toolbox.timeSteppers.',problem.timeStepper]             ;
    ts      = IntrepidTwilight.executive.buildClosure(closure,problem)  ;
    
    %   Build time-stepper closure
    closure = ['toolbox.timeSteppers.',problem.timeStepper]             ;
    ts      = IntrepidTwilight.executive.buildClosure(closure,problem)  ;


    r       = @(q) (q - ts.qLast()) - ts.deltaq(q)                      ;
    
end