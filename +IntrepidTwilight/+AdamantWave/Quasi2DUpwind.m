function q2Dup = Quasi2DUpwind(config)
    
    %   Inherit
    q2Dup = IntrepidTwilight.AdamantWave.SpaceDiscretization();
    q2Dup = q2Dup.changeID(q2Dup,'Quasi2DUpWind');


    % =================================================== %
    %                   Public Methods                    %
    % =================================================== %
    %
    q2Dup.rhs           = @(q)    rhs(q)            ;
    q2Dup.rhsMass       = @(rho)  rhsMass (rho)     ;
    q2Dup.rhsEnergy     = @(rhoe) rhsEnergy(rhoe)   ;
    q2Dup.rhsMassEnergy = @(qCV)  rhsMassEnergy(qCV);
    q2Dup.rhsMomentum   = @(rhov) rhsMomentum(rhov) ;
    %
    q2Dup.blockDiagonalJacobian = @(q) jacobianBlockDiagonal(q) ;
    %   
    q2Dup.setMass     = @(mass)     setMass(mass)           ;
    q2Dup.setMomentum = @(momentum) setMomentum(momentum)   ;
    q2Dup.setEnergy   = @(energy)   setEnergy(energy)       ;
    %
    q2Dup.update                   = @(time) update(time)           ;
    q2Dup.updateClosureEnvironment = @() updateClosureEnvironment() ;
    q2Dup.updateThermodynamicState = @() updateThermodynamicState() ;
    q2Dup.updateVelocity           = @() updateVelocity()           ;




    % =================================================== %
    %                   Self-Parameters                   %
    % =================================================== %
    q2Dup.set('epsilon'                 ,1E-8   );
    q2Dup.set('dimensionalizer.mass'    ,1      );
    q2Dup.set('dimensionalizer.energy'  ,1      );
    q2Dup.set('dimensionalizer.momentum',1      );
    if (nargin >= 1)
        q2Dup.set(config);
    end



    % =================================================== %
    %                     Model Binder                    %
    % =================================================== %
    model = [];
    q2Dup.bind = @(m) bind(m);
    function [] = bind(object)
        if isstruct(object) && object.is('model')
            model = object;
        end
    end



    % =================================================== %
    %               Declare Closure Variables             %
    % =================================================== %
    
    % ------------------------------- %
    %         Self-specified          %
    % ------------------------------- %
    
    massDim     = [];
    energyDim   = [];
    momentumDim = [];
    %
    epsilon = 1E-8;


    % ------------------------------- %
    %         Model-specified         %
    % ------------------------------- %

    mass     = [];
    energy   = [];
    momentum = [];
    %
    volCV    = [];
    %
    from  = [] ;
    to    = [] ; 
    up    = [] ;
    down  = [] ;
    %
    nCV = [] ;
    nMC = [] ;
    iM  = [] ;
    iE  = [] ;
    iP  = [] ;
    %
    volMCfrom = [] ;
    volMCto   = [] ;
    %
    Ainter = [] ;
    %
    sMass     = [] ;
    sEnergy   = [] ;
    sMomentum = [] ;
    %
    friction = [] ;
    LoD      = [] ;
    %
    volMC     = [] ;
    volMCfrom = [] ;
    volMCto   = [] ;
    %
    upDotN   =  [] ;
    downDotN =  [] ;
    %
    theta = [] ;
    g     = [] ;
    %
    Ccv    = [];
    Cmc    = [];
    Cinter = [];
    iInter = [];
    %
    massBar = [];
    vCV     = [];
    vMC     = [];
    t       = [];
    dfdq    = [];
    %
    TD = struct();
    
    
    


    % =================================================== %
    %                  Method Defintions                  %
    % =================================================== %
    q2Dup.prepare = @() prepare();
    function [] = prepare()
        assignModelValues();
        assignSelfValues();
    end
    
    function [] = assignModelValues()

        % Conserved quantities
        mass     = model.controlVolume.mass     ;
        energy   = model.controlVolume.energy   ;
        volCV    = model.controlVolume.volume   ;
        momentum = model.momentumCell.momentum  ;
        
        % Control volume sense
        from  = model.momentumCell.from ;
        to    = model.momentumCell.to   ;
        up    = model.interface.up      ;
        down  = model.interface.down    ;
        
        % Volumes of momentum cells
        volMCfrom = model.momentumCell.volumeFrom   ;
        volMCto   = model.momentumCell.volumeTo     ;
        
        % Interface parameters
        Ainter = model.interface.flowArea   ;
        
        % Sources
        sMass     = model.controlVolume.source.mass     ;
        sEnergy   = model.controlVolume.source.energy   ;
        sMomentum = model.momentumCell.source.momentum  ;
        
        % Momentum
        friction = model.momentumCell.source.friction   ;
        LoD      = model.momentumCell.LoD               ;




        % ======================================================================= %
        %                            Parameter calculating                        %
        % ======================================================================= %
        
        %   Indices
        nCV = max([from(:);to(:)])  ;
        nMC = length(from)          ;
        iM  = 1:nCV                 ;
        iE  = iM + nCV              ;
        iP  = iE + nCV              ;
        
        % Normalized (fractional volume) of momentum cells
        volMC     = volMCfrom .* volCV(from) +  volMCto .* volCV(from)  ;
        volMCfrom = volMCfrom ./ volMC                                  ;
        volMCto   = volMCto   ./ volMC                                  ;
        
        
        % Momentum cell-Interface dots
        upDotN   =  model.momentumCell.directionX(up)  .*model.interface.normalX    + ...
            model.momentumCell.directionY(up)  .*model.interface.normalY    ;
        downDotN =  model.momentumCell.directionX(down).*model.interface.normalX    + ...
            model.momentumCell.directionY(down).*model.interface.normalY    ;
        
        
        % Gravity
        theta = atan(model.momentumCell.directionY./model.momentumCell.directionX);
        g     = 9.81*cos(theta + pi/2);
        
        
        % Summation matrices
        [Ccv,Cmc,Cinter,iInter] = ...
            IntrepidTwilight.AdamantWave.toolbox.GetSummationMatrices([from,to],[up,down],[upDotN,downDotN]);
        
        
        %   Jacobi finite difference epsilon
        epsilon = 1E-8;
        
        
        % Initialization for inclusion into the closure environment
        massBar = 0;
        vCV     = 0;
        vMC     = 0;
        t       = 0;
        dfdq    = {zeros(nCV) ; zeros(nCV) ; zeros(nMC)};
        
        
        %   Create a Thermodynamic struct TD (used for passing
        %   already-calculated properties to constituitive relations)
        TD.rho  = mass   ./ volCV           ;
        TD.e    = energy ./ mass            ;
        TD.T    = Temperature(TD.rho,TD.e)  ;
        TD.P    = Pressure(TD.rho,TD.T)     ;
        TD.rhoh = TD.rho.* TD.e + TD.P      ;

    end
    
    function [] = assignSelfValues()
        epsilon     = q2Dup.get('epsilon')                  ;
        massDim     = q2Dup.get('dimensionalizer.mass')      ;
        energyDim   = q2Dup.get('dimensionalizer.energy')    ;
        momentumDim = q2Dup.get('dimensionalizer.momentum')  ;
    end

    
    
    % =================================================== %
    %                  Getters/Setters                    %
    % =================================================== %
    function [] = setMass(massStar)
        mass = massStar;
    end
    function [] = setMomentum(momentumStar)
        momentum = momentumStar;
    end
    function [] = setEnergy(energyStar)
        energy = energyStar;
    end
    
    

    % =================================================== %
    %                     Full RHS                        %
    % =================================================== %
    function f = rhs(q)
        
        % Pull conserved values
        mass     = q(iM) * massDim      ;
        energy   = q(iE) * energyDim    ;
        momentum = q(iP) * momentumDim  ;


        updateClosureEnvironment();
        
        
        f = [rhsMass()/massDim;rhsEnergy()/energyDim;rhsMomentum()/momentumDim];
        
    end


    % =================================================== %
    %             Combination Mass-Energy RHS             %
    % =================================================== %
    function f = rhsMassEnergy(q)
        
        % Pull conserved values
        mass   = q(iM) * massDim    ;
        energy = q(iE) * energyDim  ;


        f = [rhsMass()/massDim;rhsEnergy()/energyDim];
        
    end



    %{
    ===========================================================
                             Mass RHS
    ===========================================================
    %}
    function f = rhsMass(massStar)
        
        if (nargin >= 1)
            mass = massStar * massDim;
        end
        
        % Advection term
        vzRho = vCV .* (  (vCV>0).*TD.rho(from)  +  (vCV<=0).*TD.rho(to)  );
        
        % Total RHS
        f = Ccv*(vzRho) + sMass(mass,energy,momentum,TD,t) ;
        
    end



    %{
    ===========================================================
                            Energy RHS
    ===========================================================
    %}
    function f = rhsEnergy(energyStar)
        
        if (nargin >= 1)
            energy = energyStar * energyDim;
        end
        
        
        % Advection term
        vzRhoh = vCV.*(  (vCV>0).*TD.rhoh(from)  + (vCV<=0).*TD.rhoh(to)  );
        
        % Total RHS
        f  = Ccv*vzRhoh + sEnergy(mass,energy,momentum,TD,t) ;
        
    end



    %{
    ===========================================================
                          Momentum RHS
    ===========================================================
    %}
    function f = rhsMomentum(momentumStar)
        
        if (nargin >= 1)
            momentum = momentumStar * momentumDim;
        end
        
        %   Intensive momentum
        rhov = momentum ./ volMC ;
        
        % Upwind/downwind momentum advection
        fup   = rhov(up)            ;
        fdown = rhov(down)          ;
        fmom  = vMC.*(fdown - fup)  ;
        
        % Pieces
        advect = (Cmc*(fmom.*Ainter) + Cinter*(TD.P(iInter).*Ainter))   ;
        buoy   = -g.*massBar                                            ;
        fric   = -friction*LoD.*abs(rhov).*vCV                          ;
        
        % Total RHS
        f = advect + buoy + fric + sMomentum(mass,energy,momentum,TD,t) ;
        
    end







    %{
    ===========================================================
                          Update functions
    ===========================================================
    %}

    function [] = update(time)
        t = time;
    end
    
    function [] = updateClosureEnvironment()
        updateVelocity();
        updateThermodynamicState();
    end
    
    function [] = updateVelocity()
        % Volume-average density and CV surface velocities
        massBar = volMCfrom .* mass(from) + volMCto .* mass(to) ;
        vCV     = momentum ./ massBar                           ;
        
        %   Momentum cell advection
        vMC = (vCV(to) .* upDotN + vCV(down) .* downDotN) / 2;
%         denom     = 1./(massBar(up)+massBar(down))              ;
%         alphaUp   = rhov(up).*denom                             ;
%         alphaDown = rhov(up).*denom                             ;
%         vMC       = alphaUp .* upDotN + alphaDown .* downDotN   ;
    end
    
    function [] = updateThermodynamicState()
        % Thermodynamic properites
        TD.rho             = mass   ./ volCV                        ;
        TD.e               = energy ./ mass                         ;
        [TD.T,twoPhiState] = Temperature(TD.rho,TD.e,TD.T)          ;
        TD.P               = Pressure(TD.rho,TD.T,false,twoPhiState);
        TD.rhoh            = energy ./ volCV + TD.P                 ;
    end





    %{
    ===========================================================
                          Jacobians
    ===========================================================
    %}

    function dfdqOut = jacobianBlockDiagonal(q)
        
        %   Pull state values
        massRef     = q(iM) * massDim       ;
        energyRef   = q(iE) * energyDim     ;
        momentumRef = q(iP) * momentumDim   ;
        
        
        %   Perturbed values
        aboveOne = abs(massRef) >= 1;
        massStep = epsilon*(massRef.*aboveOne + not(aboveOne));
        massTD   = massRef +  massStep;
        %
        aboveOne   = abs(energyRef) >= 1;
        energyStep = epsilon*(energyRef.*aboveOne + not(aboveOne));
        energyTD   = energyRef +  energyStep;
        %
        aboveOne     = abs(momentumRef) >= 1;
        momentumStep = epsilon*(momentumRef.*aboveOne + not(aboveOne));
        momentumTD   = momentumRef +  momentumStep;
        %
        %   Build large perturbed arrays for fast Thermodynamic update
        mass     = [massRef     ; massTD      ; massRef     ; massRef       ];
        energy   = [energyRef   ; energyRef   ; energyTD    ; energyRef     ];
        momentum = [momentumRef ; momentumRef ; momentumRef ; momentumTD    ];
        
        
        %   Update and store perturbed values
        TD.T = [TD.T;TD.T;TD.T;TD.T];
        updateThermodynamicState();
        rhoBlock   = TD.rho ;
        eBlock     = TD.e   ;
        Tblock     = TD.T   ;
        Pblock     = TD.P   ;
        rhohBlock  = TD.rhoh;
        
        %   Create unmutated TD struct
        TDref.rho  = rhoBlock(1:nCV)    ;
        TDref.e    = eBlock(1:nCV)      ;
        TDref.T    = Tblock(1:nCV)      ;
        TDref.P    = Pblock(1:nCV)      ;
        TDref.rhoh = rhohBlock(1:nCV)   ;
        
        %   Store compressed perturbed vectors and reset closure variables
        massTD     = [ massRef     ; massTD     ];
        energyTD   = [ energyRef   ; energyTD   ];
        momentumTD = [ momentumRef ; momentumTD ];
        mass       = massRef                    ;
        energy     = energyRef                  ;
        momentum   = momentumRef                ;
        TD         = TDref                      ;




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
            mass(k)    = massTD(nCVk)   ;
            TD.rho(k)  = rhoBlock(nCVk) ;
            TD.e(k)    = eBlock(nCVk)   ;
            TD.T(k)    = Tblock(nCVk)   ;
            TD.P(k)    = Pblock(nCVk)   ;
            TD.rhoh(k) = rhohBlock(nCVk);
            updateVelocity();
            massPlus  = rhsMass();

            %   Calculate
            dfdq{1}(:,k) = (massPlus - mass0)/(massStep(k));
            
            %   Reset
            mass(k)    = massRef(k)     ;
            TD.rho(k)  = TDref.rho(k)   ;
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
            nCVk       = 2*nCV + k      ;
            energy(k)  = energyTD(nCV+k);
            TD.rho(k)  = rhoBlock(nCVk) ;
            TD.e(k)    = eBlock(nCVk)   ;
            TD.T(k)    = Tblock(nCVk)   ;
            TD.P(k)    = Pblock(nCVk)   ;
            TD.rhoh(k) = rhohBlock(nCVk);
            updateVelocity();
            energyPlus = rhsEnergy();

            %   Calculate
            dfdq{2}(:,k) = (energyPlus - energy0)/(energyStep(k));

            %   Reset
            energy(k)  = energyRef(k)   ;
            TD.rho(k)  = TDref.rho(k)   ;
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
            nMCk        = 3*nCV+k;
            momentum(k) = momentumTD(nMC+k) ;
            TD.rho(k)   = rhoBlock(nMCk)    ;
            TD.e(k)     = eBlock(nMCk)      ;
            TD.T(k)     = Tblock(nMCk)      ;
            TD.P(k)     = Pblock(nMCk)      ;
            TD.rhoh(k)  = rhohBlock(nMCk)   ;
            updateVelocity();
            momentumPlus = rhsMomentum();

            %   Calculate
            dfdq{3}(:,k) = (momentumPlus - momentum0)/(momentumStep(k));

            %   Reset
            momentum(k) = momentumRef(k);
            TD.rho(k)   = TDref.rho(k)  ;
            TD.e(k)     = TDref.e(k)    ;
            TD.T(k)     = TDref.T(k)    ;
            TD.P(k)     = TDref.P(k)    ;
            TD.rhoh(k)  = TDref.rhoh(k) ;

        end
        
        
        %   Reset all mutations to closure environment
        mass     = massRef      ;
        energy   = energyRef    ;
        momentum = momentumRef  ;
        TD       = TDref        ;
        
        
        %   Pass out block diagonal
        dfdqOut = dfdq;
        
    end


end