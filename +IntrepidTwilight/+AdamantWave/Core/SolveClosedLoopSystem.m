function qNDnew = SolveClosedLoopSystem(q0,s)
    
    
    q0ND = q0 ./ s.qstar;
    
    % Conserved quantities
    rho  = q0(s.rhoMask);
    rhoe = q0(s.rhoeMask);
    rhov = q0(s.rhovMask);
    Ones = ones(size(rho));
    
    DensityGuard = @(rhoe) rhoe > [996.8*Ones;Inf*Ones];
    
    fMC = @(rhov,s) s.dt*SemidiscreteUpwindMomentum(rhov,s)     ;
    fCV = @(q,s)    s.dt*SemidiscreteUpwindControlVolume(q,s)   ;
    
    % Get new momentum
    r = @(rhovNew) (rhovNew - rhov) - fMC(rhovNew,s);
    rhov = JFNKHouseholder(rhov,r,1E-8);
    
    % Get new density/energy
    r    = @(rhoeNew) (rhoeNew - [rho;rhoe]) - fCV([rhoeNew;rhov],s);
    rhoe = JFNKHouseholder([rho;rhoe],r,1E-8,DensityGuard);
    
    
    qNDnew = 0;
end

