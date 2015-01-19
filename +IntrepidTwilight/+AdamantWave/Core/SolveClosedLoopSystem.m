function q = SolveClosedLoopSystem(q,params)
    
    
    % Conserved quantities
    rhoNow   = q(1:params.nCV);
    rhoeNow  = q((params.nCV+1):(2*params.nCV));
    rhovzNow = q((2*params.nCV+1):end);
    
    
    while true
    
    % Thermodynamic properites
    params.e = rhoeNow ./ rhoNow                    ;
    params.T = Temperature(rhoNow,params.e,params.T);
    params.P = Pressure   (rhoNow,params.T)         ;
    
    
    % Update momentum field
    f = @(rhovz) SemidiscreteUpwindMomentum(rhovz,params);
    r = @(rhovz) (rhovz - rhovzNow) - params.dt*f(rhovz);
    rhovzThen = JFNKHouseholder(1.00001*rhovzNow,r,1E-8);
    
    % Pull often used values
    from  = params.from  ;
    to    = params.to    ;
    
    % Get velocities and average density
    rhoBar = params.volFrom.*rhoNow(from) + params.volTo.*rhoNow(to);
    vz     = rhovzThen ./ rhoBar ;
    rhoh   = rhoeNow + params.P ;
    
    %   Evolve control volumes
    frho     = params.Ccv*(vz.*((vz>0).* rhoNow(from) + (vz<=0).*rhoNow(to) )) + params.Srho ;
    frhoe    = params.Ccv*(vz.*((vz>0).* rhoh(from)   + (vz<=0).*rhoh(to))) + params.Srhoe;
    rhoThen  = rhoNow  + params.dt*frho  ;
    rhoeThen = rhoeNow + params.dt*frhoe ;
    
    Show(max([norm(params.dt*frho)/norm(rhoNow);norm(params.dt*frhoe)/norm(rhoeNow)]));
    
    
    rhoLimit = 0.98*MaximumDensity();
    rhoNow   = rhoThen .* (rhoThen < rhoLimit) + (rhoThen >= rhoLimit)*rhoLimit;
    rhoeNow  = rhoeThen;
    
    end
    
end
