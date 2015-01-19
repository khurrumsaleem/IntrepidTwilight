function R = Residual(q,params)
    
    rhs = SemidiscreteUpwind(q,params);
    R = (q - qOld) - dt*rhs;
    
    
end