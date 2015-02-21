function ie = ImplicitEuler(parameters)
    
    f = parameters.f;
    
    ie = @(q,qOld) f(q);
    
    
end