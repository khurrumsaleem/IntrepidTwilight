function [fValue,s] = f(q,s)
    
    % Conserved quantities
    rho   = q(1:s.Ncv);
    rhoe  = q((s.Ncv+1):(2*s.Ncv));
    rhovz = q((2*s.Ncv+1):end);
    
    % Pull often used values
    from  = s.from  ;
    to    = s.to    ;
    up    = s.up    ;
    down  = s.down  ;
    
    % Get velocities and average density
    rhoBar = s.volFrom.*rho(from) - s.volTo.*rho(to);
    vz     = rhovz ./ rhoBar ;
    vx     = vz*s.CosTheta     ;
    vy     = vz*s.SinTheta     ;
    
    % Thermodynamic properites
    e    = rhoe ./ rho                  ;
    T    = Temperature(rho,e,s.Tguess)  ;
    P    = Pressure(rho,T)              ;
    rhoh = rhoe + P                     ;
    
    % Update control volume rhs
    frho  = Ccv*(vz.*((vz>0).* rho(from)  + (vz<=0).*rho(to) )) + s.Srho ;
    frhoe = Ccv*(vz.*((vz>0).* rhoh(from) + (vz<=0).*rhoh(to))) + s.Srhoe;
    
    % Upwind/downwind momentum advection
    fup    = -rhovz(up)  .*(vx(up)  .*s.nx+vy(up)  .*s.ny)  ;
    fdown  = -rhovz(down).*(vx(down).*s.nx+vy(down).*s.ny)  ;
    fdonor = (vz>0).*fup + (vz<=0).*fdown                   ;
    
    % Update momentum rhs
    advect = -(Cinter*fdonor + Cp*P).*Ainter;
    buoy   = g*rhoBar.*cos(s.theta+pi/2)    ;
    fric   = -s.fFric*s.LoD.*abs(rhovz).*vz ;
    frhov  = advect./volume + buoy + fric   ;
    
    % Output
    fValue   = [frho;frhoe;frhov];
    s.Tguess = Tguess;

    
end