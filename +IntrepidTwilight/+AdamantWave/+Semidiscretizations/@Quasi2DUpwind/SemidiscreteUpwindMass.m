function [frho,varargout] = SemidiscreteUpwindMass(q,s)

    % Conserved quantities
    rho   = q(s.rhoMask);
    rhovz = q(s.rhovMask);
    
    % Pull often used values
    from  = s.from  ;
    to    = s.to    ;
    
    % Get velocities and average density
    rhoBar = s.volFrom.*rho(from) + s.volTo.*rho(to);
    vz     = rhovz ./ rhoBar ;
    
    % Update control volume rhs
    frho  = s.Ccv*(vz.*((vz>0).* rho(from)  + (vz<=0).*rho(to) )) + s.Srho ;
    
    if (nargout > 1)
        varargout{1} = s;
    end
    
end