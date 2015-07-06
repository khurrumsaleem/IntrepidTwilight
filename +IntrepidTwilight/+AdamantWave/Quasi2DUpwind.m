function q2Dup = Quasi2DUpwind(model)
    
    % Return closure
    q2Dup.is                    = @(s) strmpci(s,'spaceDiscretization') ;
    q2Dup.rhs                   = @(q)    rhs(q)                        ;
    q2Dup.rhsMass               = @(rho)  rhsMass (rho)                 ;
    q2Dup.rhsEnergy             = @(rhoe) rhsEnergy(rhoe)               ;
    q2Dup.rhsMomentum           = @(rhov) rhsMomentum(rhov)             ;
    q2Dup.blockDiagonalJacobian = @(q) blockDiagonalJacobian(q)         ;
    q2Dup.update                = @(time) update(time)                  ;
    
    
    
    % ======================================================================= %
    %                            Parameter unpacking                          %
    % ======================================================================= %
    
    % Conserved quantities
    rho  = model.initialState.rho0   ;
    rhoe = model.initialState.rhoe0  ;
    rhov = model.initialState.rhov0  ;
    
    rhoDim  = model.dimensionalizer.rho    ;
    rhoeDim = model.dimensionalizer.rhoe   ;
    rhovDim = model.dimensionalizer.rhov   ;
    
    % Indices
    iRho  = model.miscellaneous.iRho  ;
    iRhov = model.miscellaneous.iRhov ;
    iRhoe = model.miscellaneous.iRhoe ;
    from  = model.geometry.from       ;
    to    = model.geometry.to         ;
    up    = model.geometry.up         ;
    down  = model.geometry.down       ;
    
    % Volumes of momentum cells
    volumeBack  = model.geometry.volumeBack       ;
    volumeFront = model.geometry.volumeFront      ;
    
    % Interface parameters
    Ainter = model.geometry.Ainter   ;
    
    % Sources
    sRho  = model.miscellaneous.sRho  ;
    sRhoe = model.miscellaneous.sRhoe ;
    sRhov = model.miscellaneous.sRhov ;
    
    % Momentum 
    friction = model.miscellaneous.friction;
    LoD      = model.geometry.LoD;

    
    % ======================================================================= %
    %                            Parameter calculating                        %
    % ======================================================================= %
    
    % Amounts of stuff
    nCV    = max([from;to])             ;
    nMC    = length(from)               ;


    % Normalized (fractional volume) of momentum cells
    volTotal     = volumeBack + volumeFront;
    volumeBack   = volumeBack  ./ volTotal;
    volumeFront  = volumeFront ./ volTotal;


    % Momentum cell-Interface dots
    upDotN   =  model.geometry.zx(up)  .*model.geometry.nx  + ...
                model.geometry.zy(up)  .*model.geometry.ny  ;
    downDotN =  model.geometry.zx(down).*model.geometry.nx  + ...
                model.geometry.zy(down).*model.geometry.ny  ;


    % Gravity
    theta = atan(model.geometry.zy./model.geometry.zx);
    g     = 9.81*cos(theta + pi/2);


    % Summation matrices
    [Ccv,Cmc,Cinter,iInter] = ...
        IntrepidTwilight.AdamantWave.toolbox.GetSummationMatrices([from,to],[up,down],[upDotN,downDotN]);    


    %   Jacobi finite difference epsilon
    epsilon  = model.miscellaneous.epsilon;


    % Initialization for inclusion into the closure environment
    rhoBar = 0;
    vCV    = 0;
    vMC    = 0;
    t      = 0;
    dfdq   = {zeros(nCV) ; zeros(nCV) ; zeros(nMC)};


    %   Create a Thermodynamic struct TD (used for passing 
    %   already-calculated properties to constituitive relations)
    TD.e    = rhoe ./ rho                       ;
    TD.T    = Temperature(rho,TD.e,model.initialState.T);
    TD.P    = Pressure(rho,TD.T)                ;
    TD.rhoh = rhoe + TD.P                       ;



    function f = rhs(q)
        
        % Pull conserved values
        rho  = q(iRho)  * rhoDim    ;
        rhoe = q(iRhoe) * rhoeDim   ;
        rhov = q(iRhov) * rhovDim   ;


        updateClosureEnvironment();
        
        
        f = [rhsMass()/rhoDim;rhsEnergy()/rhoeDim;rhsMomentum()/rhovDim];
        
    end
    
    
    function [] = update(time)
        t = time;
    end
    
    function [] = updateClosureEnvironment()
        updateVelocity();
        updateThermodynamicState();
    end
    
    function [] = updateVelocity()
        % Get average density and CV surface velocities
        rhoBar    = volumeBack.*rho(from) + volumeFront.*rho(to);
        
        %   Control volume advection
        vCV       = rhov ./ rhoBar;
        
        %   Momentum cell advection
        denom     = 1./(rhoBar(up)+rhoBar(down))                ;
        alphaUp   = rhov(up).*denom                             ;
        alphaDown = rhov(up).*denom                             ;
        vMC       = alphaUp .* upDotN + alphaDown .* downDotN   ;
    end
    
    function [] = updateThermodynamicState()
        % Thermodynamic properites
        TD.e               = rhoe ./ rho                            ;
        [TD.T,twoPhiState] = Temperature(rho,TD.e,TD.T)             ;
        TD.P               = Pressure(rho,TD.T,false,twoPhiState)   ;
        TD.rhoh            = rhoe + TD.P                            ;
    end






    %{
    ===========================================================
                             Mass RHS
    ===========================================================
    %}
    function f = rhsMass(rhoStar)
        
        if (nargin >= 1)
            rho = rhoStar * rhoDim;
        end
        
        % Advection term
        vzRho = vCV.*(  (vCV>0).*rho(from) + (vCV<=0).*rho(to)  );
        
        % Total RHS
        f = Ccv*(vzRho) + sRho(rho,rhoe,rhov,TD,t) ;
        
    end

    
    
    
    %{
    ===========================================================
                            Energy RHS
    ===========================================================
    %}
    function f = rhsEnergy(rhoeStar)
        
        if (nargin >= 1)
            rhoe = rhoeStar * rhoeDim;
        end
        
        
        % Advection term
        vzRhoh = vCV.*(  (vCV>0).*TD.rhoh(from)  + (vCV<=0).*TD.rhoh(to)  );
        
        % Total RHS
        f  = Ccv*vzRhoh + sRhoe(rho,rhoe,rhov,TD,t) ;
        
    end
    
    

    %{
    ===========================================================
                          Momentum RHS
    ===========================================================
    %}
    function f = rhsMomentum(rhovStar)
        
        if (nargin >= 1)
            rhov = rhovStar * rhovDim;
        end
        
        
        % Upwind/downwind momentum advection
        fup   = rhov(up)            ;
        fdown = rhov(down)          ;
        fmom  = vMC.*(fdown - fup)  ;
        
        % Pieces
        advect = (Cmc*(fmom.*Ainter) + Cinter*(TD.P(iInter).*Ainter));
        buoy   = -g.*rhoBar                     ;
        fric   = -friction*LoD.*abs(rhov).*vCV  ;
        
        % Total RHS
        f = advect + buoy + fric + sRhov(rho,rhoe,rhov,TD,t) ;
        
    end
    
    
    
    function dfdqOut = blockDiagonalJacobian(q)
        
        %   Pull state values
        rhoRef  = q(iRho)  * rhoDim     ;
        rhoeRef = q(iRhoe) * rhoeDim    ;
        rhovRef = q(iRhov) * rhovDim    ;
        
        
        %   Perturbed values
        aboveOne = abs(rhoRef) >= 1;
        rhoStep  = epsilon*(rhoRef.*aboveOne + not(aboveOne));
        rhoTD    = rhoRef +  rhoStep;
        %
        aboveOne = abs(rhoeRef) >= 1;
        rhoeStep  = epsilon*(rhoeRef.*aboveOne + not(aboveOne));
        rhoeTD    = rhoeRef +  rhoeStep;
        %
        aboveOne = abs(rhovRef) >= 1;
        rhovStep  = epsilon*(rhovRef.*aboveOne + not(aboveOne));
        rhovTD    = rhovRef +  rhovStep;
        %
        %   Build large perturbed arrays for fast Thermodynamic update
        rho    = [rhoRef  ; rhoTD   ; rhoRef  ; rhoRef  ];
        rhoe   = [rhoeRef ; rhoeRef ; rhoeTD  ; rhoeRef ];
        rhov   = [rhovRef ; rhovRef ; rhovRef ; rhovTD  ];
        
        
        %   Update and store perturbed values
        TD.T = [TD.T;TD.T;TD.T;TD.T];
        updateThermodynamicState();
        Tblock     = TD.T    ;
        Pblock     = TD.P    ;
        rhohBlock  = TD.rhoh ;
        
        %   Create unmutated TD struct
        TDref.e    = rhoeRef ./ rhoRef  ;
        TDref.T    = Tblock(1:nCV)      ;
        TDref.P    = Pblock(1:nCV)      ;
        TDref.rhoh = rhohBlock(1:nCV)   ;
        
        %   Store compressed perturbed vectors and reset closure variables
        rhoTD  = [ rhoRef  ; rhoTD  ]   ;
        rhoeTD = [ rhoeRef ; rhoeTD ]   ;
        rhovTD = [ rhovRef ; rhovTD ]   ;
        rho    = rhoRef                 ;
        rhoe   = rhoeRef                ;
        rhov   = rhovRef                ;
        TD     = TDref                  ;




        % ========================================================= %
        %                         Mass Block                        %
        % ========================================================= %
        
        %   Get unperturbed value
        updateVelocity();
        mass0 = rhsMass();


        % Iterate through columns of the block
        K = 1:nCV;
        for k = K

            %   Perturb
            nCVk       = nCV+k;
            rho(k)     = rhoTD(nCVk);
            TD.e(k)    = rhoeRef(k) / rho(k);
            TD.T(k)    = Tblock(nCVk);
            TD.P(k)    = Pblock(nCVk);
            TD.rhoh(k) = rhohBlock(nCVk);
            updateVelocity();
            massPlus  = rhsMass();

            %   Calculate
            dfdq{1}(:,k) = (massPlus - mass0)/(rhoStep(k));
            
            %   Reset
            rho(k)     = rhoRef(k)      ;
            TD.e(k)    = TDref.e(k)     ;
            TD.T(k)    = TDref.T(k)     ;
            TD.P(k)    = TDref.P(k)     ;
            TD.rhoh(k) = TDref.rhoh(k)  ;

        end



        % ========================================================= %
        %                       Energy Block                        %
        % ========================================================= %
        
        %   Get unperturbed value
        updateVelocity();
        energy0 = rhsEnergy();
        
        %   Iterate through columns of the block
        for k = K

            %   Perturb
            nCVk       = 2*nCV+k;
            rhoe(k)    = rhoeTD(nCV+k);
            TD.e(k)    = rhoe(k) / rho(k);
            TD.T(k)    = Tblock(nCVk);
            TD.P(k)    = Pblock(nCVk);
            TD.rhoh(k) = rhohBlock(nCVk);
            updateVelocity();
            energyPlus = rhsEnergy();

            %   Calculate
            dfdq{2}(:,k) = (energyPlus - energy0)/(rhoeStep(k));

            %   Reset
            rhoe(k)    = rhoeRef(k)     ;
            TD.e(k)    = TDref.e(k)     ;
            TD.T(k)    = TDref.T(k)     ;
            TD.P(k)    = TDref.P(k)     ;
            TD.rhoh(k) = TDref.rhoh(k)  ;

        end



        % ========================================================= %
        %                       Momentum Block                      %
        % ========================================================= %
        
        %   Get unperturbed value
        updateVelocity();
        momentum0 = rhsMomentum();

        %   Iterate through columns of the block
        K = 1:nMC;
        for k = K
            
            %   Perturb
            nMCk       = 3*nCV+k;
            rhov(k)    = rhovTD(nMC+k);
            TD.T(k)    = Tblock(nMCk);
            TD.P(k)    = Pblock(nMCk);
            TD.rhoh(k) = rhohBlock(nMCk);
            updateVelocity();
            momentumPlus = rhsMomentum();

            %   Calculate
            dfdq{3}(:,k) = (momentumPlus - momentum0)/(rhovStep(k));

            %   Reset
            rhov(k)    = rhovRef(k)     ;
            TD.e(k)    = TDref.e(k)     ;
            TD.T(k)    = TDref.T(k)     ;
            TD.P(k)    = TDref.P(k)     ;
            TD.rhoh(k) = TDref.rhoh(k)  ;

        end
        
        
        %   Reset all mutations to closure environment
        rho  = rhoRef   ;
        rhoe = rhoeRef  ;
        rhov = rhovRef  ;
        TD   = TDref    ;
        
        
        %   Pass out block diagonal
        dfdqOut = dfdq;
        
    end
    
    
    
    
    
    
end