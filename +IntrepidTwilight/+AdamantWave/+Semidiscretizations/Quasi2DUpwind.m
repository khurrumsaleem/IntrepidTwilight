function q2Dup = Quasi2DUpwind(problem)
    
    % Return closure
    q2Dup.rhs                   = @(q) rhs(q);
    q2Dup.blockDiagonalJacobian = @(q) blockDiagonalJacobian(q);
    
    
    
    % ======================================================================= %
    %                            Parameter unpacking                          %
    % ======================================================================= %
    
    % Conserved quantities
    rho  = problem.initialState.rho0   ;
    rhoe = problem.initialState.rhoe0  ;
    rhov = problem.initialState.rhov0  ;
    
    rhoDim  = problem.dimensionalizer.rho    ;
    rhoeDim = problem.dimensionalizer.rhoe   ;
    rhovDim = problem.dimensionalizer.rhov   ;
    
    % Indices
    iRho  = problem.miscellaneous.iRho  ;
    iRhov = problem.miscellaneous.iRhov ;
    iRhoe = problem.miscellaneous.iRhoe ;
    from  = problem.geometry.from       ;
    to    = problem.geometry.to         ;
    up    = problem.geometry.up         ;
    down  = problem.geometry.down       ;
    
    % Volumes of momentum cells
    volumeBack  = problem.geometry.volumeBack       ;
    volumeFront = problem.geometry.volumeFront      ;
    
    % Interface parameters
    Ainter = problem.geometry.Ainter   ;
    
    % Sources
    sRho  = problem.miscellaneous.sRho  ;
    sRhoe = problem.miscellaneous.sRhoe ;
    sRhov = problem.miscellaneous.sRhov ;
    
    % Momentum 
    friction = problem.miscellaneous.friction;
    LoD      = problem.geometry.LoD;

    
    % ======================================================================= %
    %                            Parameter calculating                        %
    % ======================================================================= %
    
    % Amounts of stuff
    nCV    = max([from;to])             ;
    nMC    = length(from)               ;
%     nInter = length(problem.geometry.nx);
%     nEq    = 2*nCV + nMC                ;


    % Normalized (fractional volume) of momentum cells
    volTotal     = volumeBack + volumeFront;
    volumeBack   = volumeBack  ./ volTotal;
    volumeFront  = volumeFront ./ volTotal;


    % Momentum cell-Interface dots
    upDotN   =  problem.geometry.zx(up)  .*problem.geometry.nx  + ...
                problem.geometry.zy(up)  .*problem.geometry.ny  ;
    downDotN =  problem.geometry.zx(down).*problem.geometry.nx  + ...
                problem.geometry.zy(down).*problem.geometry.ny  ;


    % Gravity
    theta = atan(problem.geometry.zy./problem.geometry.zx);
    g     = 9.81*cos(theta + pi/2);


    % Summation matrices
    [Ccv,Cmc,Cinter,iInter] = ...
        IntrepidTwilight.AdamantWave.toolbox.GetSummationMatrices([from,to],[up,down],[upDotN,downDotN]);    


    %   Jacobi finite difference epsilon
    epsilon  = problem.miscellaneous.epsilon;
    plusEps  = 1 + epsilon;
    minusEps = 1 - epsilon;


    % Initialization for inclusion into the closure environment
    rhoBar = 0;
    v      = 0;


    %   Create a Thermodynamic struct TD (used for passing 
    %   already-calculated properties to constituitive relations)
    TD.e    = rhoe ./ rho                       ;
    TD.T    = Temperature(rho,TD.e,problem.initialState.T);
    TD.P    = Pressure(rho,TD.T)                ;
    TD.rhoh = rhoe + TD.P                       ;



    function f = rhs(q)
        
        % Pull conserved values
        rho  = q(iRho)  * rhoDim    ;
        rhoe = q(iRhoe) * rhoeDim   ;
        rhov = q(iRhov) * rhovDim   ;


        updateClosureEnvironment();
        
        
        f = [massRHS()/rhoDim;energyRHS()/rhoeDim;momentumRHS()/rhovDim];
        
    end
    
    
    function [] = updateClosureEnvironment()
        updateVelocity();
        updateThermodynamicState();
    end
    
    function [] = updateVelocity()
        % Get average density and CV surface velocities
        rhoBar = volumeBack.*rho(from) + volumeFront.*rho(to);
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
        rhoRef  = q(iRho)  * rhoDim     ;
        rhoeRef = q(iRhoe) * rhoeDim    ;
        rhovRef = q(iRhov) * rhovDim    ;
        
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
            dfdqBD(K,k) = (massPlus - massMinus)/(2*epsilon*rhoDim);

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
            dfdqBD(nCV+K,k) = (energyPlus - energyMinus)/(2*epsilon*rhoeDim);


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
            dfdqBD(2*nCV+K,k) = (momentumPlus - momentumMinus)/(2*epsilon*rhovDim);

        end
        
        
        % Reset all mutations to closue enviroment
        rho  = rhoRef  / rhoDim     ;
        rhoe = rhoeRef / rhoeDim    ;
        rhov = rhovRef / rhovDim    ;
        TD   = TDref    ;
        
    end
    
    
end