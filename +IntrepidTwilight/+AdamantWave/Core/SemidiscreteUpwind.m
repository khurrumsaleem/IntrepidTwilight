function f = SemidiscreteUpwind(q,s)
    
    % Conserved quantities
    rho   = q(s.rhoMask);
    rhoe  = q(s.rhoeMask);
    rhovz = q(s.rhovMask);
    
    % Pull often used values
    from  = s.from  ;
    to    = s.to    ;
    up    = s.up    ;
    down  = s.down  ;
    
    % Get velocities and average density
    rhoBar = s.volFrom.*rho(from) + s.volTo.*rho(to);
    vz     = rhovz ./ rhoBar ;
    vavg   = (vz(up).*s.upDotN + vz(down).*s.downDotN)/2;
    
    % Thermodynamic properites
    e    = rhoe ./ rho              ;
    T    = Temperature(rho,e);
    P    = Pressure(rho,T)          ;
    rhoh = rhoe + P                 ;
    
    % Update control volume rhs
    frho  = s.Ccv*(vz.*((vz>0).* rho(from)  + (vz<=0).*rho(to) )) + s.Srho ;
    frhoe = s.Ccv*(vz.*((vz>0).* rhoh(from) + (vz<=0).*rhoh(to))) + s.Srhoe;
    
    % Upwind/downwind momentum advection
    fup    = rhovz(up)  .*vz(up)  .*s.upDotN    ;
    fdown  = rhovz(down).*vz(down).*s.downDotN  ;
    fdonor = (vavg>0).*fup + (vavg<=0).*fdown   ;
    
    % Update momentum rhs
    advect = (s.Cmc*(fdonor.*s.Ainter) + s.Cp*(P(s.inter).*s.Ainter));
    buoy   = s.g.*rhoBar                    ;
    fric   = -s.fFric*s.LoD.*abs(rhovz).*vz ;
    frhov  = advect + buoy + fric           ;
    
    % Output
    f = [frho;frhoe;frhov];
    f = f(s.internalMask);

    
end