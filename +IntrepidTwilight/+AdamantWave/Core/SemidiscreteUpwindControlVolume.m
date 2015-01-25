function fCV = SemidiscreteUpwindControlVolume(q,s)
    
    % Conserved quantities
    rho   = q(s.rhoMask);
    rhoe  = q(s.rhoeMask);
    rhovz = q(s.rhovMask);
    
    % Pull often used values
    from  = s.from  ;
    to    = s.to    ;
    
    % Get velocities and average density
    rhoBar = s.volFrom.*rho(from) + s.volTo.*rho(to);
    vz     = rhovz ./ rhoBar ;
    
    % Thermodynamic properites
    e    = rhoe ./ rho          ;
    T    = Temperature(rho,e)   ;
    P    = Pressure(rho,T)      ;
    rhoh = rhoe + P             ;
    
    % Update control volume rhs
    frho  = s.Ccv*(vz.*((vz>0).* rho(from)  + (vz<=0).*rho(to) )) + s.Srho ;
    frhoe = s.Ccv*(vz.*((vz>0).* rhoh(from) + (vz<=0).*rhoh(to))) + s.Srhoe;
    
    % Output
    fCV = [frho;frhoe];

    
end