function q2Dup = Quasi2DUpwind(parameters)
    
    % Return closure
    q2Dup.flux = @(q) flux(q);
    
    
    % Conserved quantities
    rho  = parameters.rho0   ;
    rhoe = parameters.rhoe0  ;
    rhov = parameters.rhov0  ;
    
    % Indices
    iRho  = parameters.iRho  ;
    iRhov = parameters.iRhov ;
    iRhoe = parameters.iRhoe ;
    from  = parameters.from  ;
    to    = parameters.to    ;
    back  = parameters.back  ;
    front = parameters.front ;
    up    = parameters.up    ;
    down  = parameters.down  ;
    
    % Dot products
    upDotN   = parameters.upDotN     ;
    downDotN = parameters.downDotN   ;
    
    % Connectivity/Summation matrices
    Ccv  = parameters.Ccv    ;
    Cmc  = parameters.Cmc    ;
    
    % Interface parameters
    Cinter = parameters.Cinter   ;
    iInter = parameters.iInter   ;
    Ainter = parameters.Ainter   ;
    
    % Sources
    sRho  = parameters.sRho  ;
    sRhoe = parameters.sRhoe ;
    sRhov = parameters.sRhov ;
    
    % Momentum 
    friction = parameters.friction;
    LoD      = parameters.LoD;
    g        = parameters.g;
    
    
    
    % Explicit unitialization for closure
    rhoBar = 0;
    v      = 0;
    
    %   Create a Thermodynamic struct TD (used for passing already
    %   calculated properties to constituitive relations)
    TD.e    = rhoe ./ rho               ;
    TD.T    = Temperature(rho,TD.e,parameters.T) ;
    TD.P    = Pressure(rho,TD.T)        ;
    TD.rhoh = rhoe + TD.P               ;
    
    
    function f = flux(q)
        
        % Pull conserved values
        rho  = q(iRho)	;
        rhoe = q(iRhoe) ;
        rhov = q(iRhov) ;
        
        
        updateClosureEnvironment();
        
        
        f = [mass();energy();momentum()];
        
    end
    
    
    function [] = updateClosureEnvironment()
        updateVelcoity();
        updateThermodynamics();
    end
    
    function [] = updateVelcoity()
        % Get average density and CV surface velocities
        rhoBar = back.*rho(from) + front.*rho(to);
        v     = rhov ./ rhoBar ;
    end
    
    function [] = updateThermodynamics()
        % Thermodynamic properites
        TD.e    = rhoe ./ rho               ;
        TD.T    = Temperature(rho,TD.e,TD.T);
        TD.P    = Pressure(rho,TD.T)        ;
        TD.rhoh = rhoe + TD.P               ;
    end






    %{
    ===========================================================
                             Mass RHS
    ===========================================================
    %}
    function f = mass()
        
        % Advection term
        vzRho = v.*(  (v>0).*rho(from) + (v<=0).*rho(to)  );
        
        % Total RHS
        f = Ccv*(vzRho) + sRho ;
        
    end
    
    
    
    
    
    %{
    ===========================================================
                            Energy RHS
    ===========================================================
    %}
    function f = energy()
        
        % Advection term
        vzRhoh = v.*(  (v>0).*TD.rhoh(from)  + (v<=0).*TD.rhoh(to)  );
        
        % Total RHS
        f  = Ccv*vzRhoh + sRhoe ;
        
    end
    
    

    %{
    ===========================================================
                          Momentum RHS
    ===========================================================
    %}
    function f = momentum()
        
        
        vavg   = (v(up).*upDotN + v(down).*downDotN)/2;
        
        % Upwind/downwind momentum advection
        fup   = rhov(up)  .*(v(up)  .*upDotN)   ;
        fdown = rhov(down).*(v(down).*downDotN) ;
        fmom  = sign(vavg).*(fdown - fup)       ;
        
        % Pieces
        advect = (Cmc*(fmom.*Ainter) + Cinter*(TD.P(iInter).*Ainter));
        buoy   = -g.*rhoBar                   ;
        fric   = -friction*LoD.*abs(rhov).*v ;
        
        % Total RHS
        f = advect + buoy + fric + sRhov ;
        
    end
    
end