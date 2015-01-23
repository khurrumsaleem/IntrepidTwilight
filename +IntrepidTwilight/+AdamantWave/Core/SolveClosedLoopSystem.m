function qNDnew = SolveClosedLoopSystem(q0,s)
    
    qND = q0 ./ s.qstar;
    
    f         = @(q) SemidiscreteUpwind(q.*s.qstar,s);
    fInternal = @(q) s.dt*f(q)./s.qstar(s.internalMask);
    rINternal = internalResidual(qND,fInternal,s.internalMask);
    rBoundary = @(q) [q(1) - q(s.nCV);q(s.nCV+1) - q(2*s.nCV)];
    r         = @(q)[rBoundary(q);rINternal(q)];
    
    qNDnew = JFNKHouseholder(1.0001*qND,r,1E-8);
    
end

function r = internalResidual(q0,f,mask)
    r = @(q) q(mask) - q0(mask) - f(q);
end