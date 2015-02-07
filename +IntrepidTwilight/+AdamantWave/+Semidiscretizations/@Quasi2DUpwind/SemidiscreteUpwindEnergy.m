function [frhoe,varargout] = SemidiscreteUpwindEnergy(q,s)

    % Conserved quantities
    rho   = q(s.rhoMask);
    rhoe  = q(s.rhoeMask);
    rhovz = q(s.rhovMask);
    
    % Pull often used values
    from  = s.from  ;
    to    = s.to    ;
    
    % Get velocities and average density
    rhoBar = s.volFrom.*rhoe(from) + s.volTo.*rhoe(to);
    vz     = rhovz ./ rhoBar ;
    
    % Thermodynamic properites
    e    = rhoe ./ rho              ;
    T    = Temperature(rho,e,s.T)   ;
    P    = Pressure(rho,T)          ;
    rhoh = rhoe + P                 ;
    
    % Update control volume rhs
    frhoe  = s.Ccv*(vz.*((vz>0).* rhoh(from)  + (vz<=0).*rhoh(to) )) + s.Srhoe ;
    
    if (nargout > 1)
        varargout{1} = s;
    end
    
end