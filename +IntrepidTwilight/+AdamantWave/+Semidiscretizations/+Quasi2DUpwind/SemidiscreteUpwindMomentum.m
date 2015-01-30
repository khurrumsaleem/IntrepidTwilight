function [frhov,varargout] = SemidiscreteUpwindMomentum(q,s)
    
    % Conserved quantities
    rho   = q(s.rhoMask);
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
    
    % Upwind/downwind momentum advection
    fup   = rhovz(up)  .*(vz(up)  .*s.upDotN)   ;
    fdown = rhovz(down).*(vz(down).*s.downDotN) ;
    fmom  = sign(vavg).*(fdown - fup)           ;
    
    % Update momentum rhs
    advect = (s.Cmc*(fmom.*s.Ainter) + s.Cp*(s.P(s.inter).*s.Ainter));
    buoy   = -s.g.*rhoBar                   ;
    fric   = -s.fFric*s.LoD.*abs(rhovz).*vz ;
    frhov  = advect + buoy + fric           ;
    
    if (nargout > 1)
        varargout{1} = s;
    end
    
end