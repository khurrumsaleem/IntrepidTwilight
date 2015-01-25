function frhov = SemidiscreteUpwindMomentum(rhovz,s)
    
    % Pull often used values
    from  = s.from  ;
    to    = s.to    ;
    up    = s.up    ;
    down  = s.down  ;
    
    % Get velocities and average density
    rhoBar = s.volFrom.*s.rho(from) + s.volTo.*s.rho(to);
    vz     = rhovz ./ rhoBar ;
    vavg   = (vz(up).*s.upDotN + vz(down).*s.downDotN)/2;
    
    % Upwind/downwind momentum advection
    fup    = rhovz(up)  .*vz(up)  .*s.upDotN    ;
    fdown  = rhovz(down).*vz(down).*s.downDotN  ;
    fdonor = (vavg>0).*fup + (vavg<=0).*fdown   ;
    
    % Update momentum rhs
    advect = (s.Cmc*(fdonor.*s.Ainter) + s.Cp*(s.P(s.inter).*s.Ainter));
    buoy   = s.g.*rhoBar                    ;
    fric   = -s.fFric*s.LoD.*abs(rhovz).*vz ;
    frhov  = advect + buoy + fric           ;
        
end