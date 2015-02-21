function q2Dup = Quasi2DUpwind(parameters)
    
    % Return closure
    q2Dup.rhs                   = @(q) rhs(q);
    q2Dup.blockDiagonalJacobian = @(q) blockDiagonalJacobian(q);
    
    
    
    nCV = parameters.nCV;
    nMC = parameters.nMC;
    
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
    
    % 
    epsilon  = parameters.epsilon;
    plusEps  = 1 + epsilon;
    minusEps = 1 - epsilon;
    
    % Initialization for inclusion into the closure environment
    rhoBar = 0;
    v      = 0;


    %   Create a Thermodynamic struct TD (used for passing 
    %   already-calculated properties to constituitive relations)
    TD.e    = rhoe ./ rho                       ;
    TD.T    = Temperature(rho,TD.e,parameters.T);
    TD.P    = Pressure(rho,TD.T)                ;
    TD.rhoh = rhoe + TD.P                       ;
    
    
    function f = rhs(q)
        
        % Pull conserved values
        rho  = q(iRho)	;
        rhoe = q(iRhoe) ;
        rhov = q(iRhov) ;
        
        
        updateClosureEnvironment();
        
        
        f = [massRHS();energyRHS();momentumRHS()];
        
    end
    
    
    function [] = updateClosureEnvironment(rho,rhoe,rhov)
        updateVelocity(rho,rhoe,rhov);
        updateThermodynamicState(rho,rhoe,rhov);
    end
    
    function [] = updateVelocity()
        % Get average density and CV surface velocities
        rhoBar = back.*rho(from) + front.*rho(to);
        v     = rhov ./ rhoBar ;
    end
    
    function [] = updateThermodynamicState()
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
    function f = massRHS()
        
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
    function f = energyRHS()
        
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
    function f = momentumRHS()
        
        
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
    
    
    
    function dfdqBD = blockDiagonalJacobian(q)
        
        %   Pull state values
        rhoRef  = q(iRho)  ;
        rhoeRef = q(iRhoe) ;
        rhovRef = q(iRhov) ;
        
        %   Build large perturbed arrays for fast Thermodynamic update
        rhoTD  = [ plusEps*rhoRef  ; minusEps*rhoRef  ];
        rhoeTD = [ plusEps*rhoeRef ; minusEps*rhoeRef ];
        rhovTD = [ plusEps*rhovRef ; minusEps*rhovRef ];
        rho    = [rhoRef ;      rhoTD      ;  rhoRef;rhoRef  ;  rhoRef;rhoRef  ];
        rhoe   = [rhoeRef; rhoeRef;rhoeRef ;     rhoeTD      ; rhoeRef;rhoeRef ];
        rhov   = [rhovRef; rhovRef;rhovRef ; rhovRef;rhovRef ;     rhovTD      ];
        
        
        %   Update and store perturbed values
        TD.T = [TD.T;TD.T;TD.T;TD.T;TD.T;TD.T;TD.T];
        updateThermodynamicState();
        Tblock     = TD.T    ;
        Pblock     = TD.P    ;
        rhohBlock  = TD.rhoh ;
        
        %   Create unmutated TD struct
        TDref.e    = rhoeRef ./ rhoRef ;
        TDref.T    = Tblock(1:nCV)     ;
        TDref.P    = Tblock(1:nCV)     ;
        TDref.rhoh = Tblock(1:nCV)     ;
        
        %   Store compressed perturbed vectors and reset closure variables
        rhoTD  = [ rhoRef  ; rhoTD  ]   ;
        rhoeTD = [ rhoeRef ; rhoeTD ]   ;
        rhovTD = [ rhovRef ; rhovTD ]   ;
        rho    = rhoRef                 ;
        rhoe   = rhoeRef                ;
        rhov   = rhovRef                ;


        dfdqBD = zeros([2*nCV+nMC,max(nCV,nMC)]);

        % ========================================================= %
        %                         Mass Block                        %
        % ========================================================= %
        K = 1:nCV;
        for k = K
            
            %   Reset rho and TD
            rho = rhoRef    ;
            TD  = TDref     ;
            
            %   Plus perturbation
            nCVk       = nCV+k;
            rho(k)     = rhoTD(nCVk);
            TD.e(k)    = rhoeRef(k) / rho(k);
            TD.T(k)    = Tblock(nCVk);
            TD.P(k)    = Pblock(nCVk);
            TD.rhoh(k) = rhohBlock(nCVk);
            updateVelocity();
            massPlus  = massRHS();
            
            %   Minus perturbation
            nCVk       = 2*nCV+k;
            rho(k)     = rhoTD(nCVk);
            TD.e(k)    = rhoeRef(k) / rhoTD(nCVk);
            TD.T(k)    = Tblock(nCVk);
            TD.P(k)    = Pblock(nCVk);
            TD.rhoh(k) = rhohBlock(nCVk);
            updateVelocity();
            massMinus = massRHS();
            
            %   Calculate
            dfdqBD(K,k) = (massPlus - massMinus)/(2*epsilon);

        end
        %   Reset rho
        rho = rhoRef;


        
        % ========================================================= %
        %                       Energy Block                        %
        % ========================================================= %
        for k = K
            
            
            %   Reset rhoe and TD
            rhoe = rhoeRef  ;
            TD   = TDref    ;
            
            %   Plus perturbation
            nCVk       = 3*nCV+k;
            rhoe(k)    = rhoeTD(nCV+k);
            TD.e(k)    = rhoe(k) / rho(k);
            TD.T(k)    = Tblock(nCVk);
            TD.P(k)    = Pblock(nCVk);
            TD.rhoh(k) = rhohBlock(nCVk);
            updateVelocity();
            energyPlus = energyRHS();
            
            %   Minus perturbation
            nCVk       = 4*nCV+k;
            rhoe(k)    = rhoeTD(2*nCV+k);
            TD.e(k)    = rhoe(k) / rho(k);
            TD.T(k)    = Tblock(nCVk);
            TD.P(k)    = Pblock(nCVk);
            TD.rhoh(k) = rhohBlock(nCVk);
            updateVelocity();
            energyMinus = energyRHS();
            
            %   Calculate
            dfdqBD(nCV+K,k) = (energyPlus - energyMinus)/(2*epsilon);


        end
        %   Reset rho
        rhoe = rhoeRef;


        % ========================================================= %
        %                       Momentum Block                      %
        % ========================================================= %
        K = 1:nMC;
        for k = K
            
            %   Reset rhov and TD
            rhov = rhovRef;
            TD   = TDref;
            
            %   Plus perturbation
%             nMCk    = 5*nCV+k;
            rhov(k)    = rhovTD(nMC+k);
%             TD.T(k)    = Tblock(nMCk);
%             TD.P(k)    = Pblock(nMCk);
%             TD.rhoh(k) = rhohBlock(nMCk);
            updateVelocity();
            momentumPlus = momentumRHS();
            
            %   Minus perturbation
%             nMCk       = 6*nCV+k;
            rhov(k)    = rhovTD(2*nMC+k);
%             TD.T(k)    = Tblock(nMCk);
%             TD.P(k)    = Pblock(nMCk);
%             TD.rhoh(k) = rhohBlock(nMCk);
            updateVelocity();
            momentumMinus = momentumRHS();
            
            
            %   Calculate
            dfdqBD(2*nCV+K,k) = (momentumPlus - momentumMinus)/(2*epsilon);

        end
        
        
        % Reset all mutations to closue enviroment
        rho  = rhoRef   ;
        rhoe = rhoeRef  ;
        rhov = rhovRef  ;
        TD   = TDref    ;
        
    end
    
    
end