function ie = ImplicitEuler(parameters)
    
    f = parameters.f;
    
    ie.deltaq = @(q) f(q);
    ie.qOldUpdate = @(qold) [] ;
    
end