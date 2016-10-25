function q2Dup = Quasi2DUpwind(config)
    
    %   Inherit
    q2Dup = IntrepidTwilight.AdamantWave.SpaceDiscretization();
    q2Dup = q2Dup.changeID(q2Dup,'Quasi2DUpWind');
    
    
    % =================================================== %
    %                   Public Methods                    %
    % =================================================== %
    %
    q2Dup.rhs           = @(q,t)    rhs(q,t)            ;
    q2Dup.rhsMass       = @(rho,t)  rhsMass (rho,t)     ;
    q2Dup.rhsEnergy     = @(rhoe,t) rhsEnergy(rhoe,t)   ;
    q2Dup.rhsVolume     = @(vol,t)  rhsVolume(vol,t)    ;
    q2Dup.rhsMomentum   = @(rhov,t) rhsMomentum(rhov,t) ;
    %
    q2Dup.rhsFast = @(q,t) rhsFast(q,t);
    q2Dup.rhsSlow = @(q,t) rhsSlow(q,t);
    %
    q2Dup.makeDimensional       = @(q) makeDimensional(q)       ;
    q2Dup.makeFastDimensional   = @(q) makeFastDimensional(q)   ;
    q2Dup.makeSlowDimensional   = @(q) makeSlowDimensional(q)   ;
    q2Dup.makeDimensionless     = @(q) makeDimensionless(q)     ;
    q2Dup.makeFastDimensionless = @(q) makeFastDimensionless(q) ;
    q2Dup.makeSlowDimensionless = @(q) makeSlowDimensionless(q) ;
    %
    q2Dup.blockDiagonalJacobian = @(q) jacobianBlockDiagonal(q) ;
    %
    q2Dup.setMass     = @(mass)     setMass(mass)           ;
    q2Dup.setMomentum = @(momentum) setMomentum(momentum)   ;
    q2Dup.setEnergy   = @(energy)   setEnergy(energy)       ;
    q2Dup.setVolume   = @(volume)   setVolume(volume)       ;
    q2Dup.setAll      = @(q) setAll(q)                      ;
    q2Dup.setFast     = @(q) setFast(q)                     ;
    q2Dup.setSlow     = @(q) setSlow(q)                     ;
    %
    q2Dup.getMass     = @() getMass()       ;
    q2Dup.getEnergy   = @() getEnergy()     ;
    q2Dup.getVolume   = @() getVolume()     ;
    q2Dup.getMomentum = @() getMomentum()   ;
    q2Dup.getAll      = @() getAll()        ;
    q2Dup.getFast     = @() getFast()       ;
    q2Dup.getSlow     = @() getSlow()       ;
    %
    q2Dup.updateClosureEnvironment = @() updateClosureEnvironment() ;
    q2Dup.updateThermodynamicState = @() updateThermodynamicState() ;
    q2Dup.updateVelocity           = @() updateVelocity()           ;
    
    
    
    
    % =================================================== %
    %                   Self-Parameters                   %
    % =================================================== %
    q2Dup.set('epsilon'                 , 1E-6  );
    q2Dup.set('dimensionalizer.mass'    , 1     );
    q2Dup.set('dimensionalizer.energy'  , 1     );
    q2Dup.set('dimensionalizer.momentum', 1     );
    q2Dup.set('isDynamicVolume'         , true  );
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
    %        Internal use only        %
    % ------------------------------- %
    isNotPrepared = true;
    
    
    % ------------------------------- %
    %         Self-specified          %
    % ------------------------------- %
    
    massDim     = [];
    energyDim   = [];
    volumeDim   = [];
    momentumDim = [];
    %
    epsilon = 1E-8;
    
    
    % ------------------------------- %
    %         Model-specified         %
    % ------------------------------- %
    
    %   Conserved quantities
    mass     = [];
    energy   = [];
    volume   = [];
    momentum = [];
    
    
    %   Velocities
    vCV = [];
    vI  = [];
    
    %   Connectivity Arrays: Control Volume sense
    from  = [] ;
    to    = [] ;
    
    %   Connectivity Arrays: Momentum Interface sense
    up    = [] ;
    down  = [] ;
    
    %   Indexing arrays
    nCV = [] ;
    nMC = [] ;
    iM  = [] ;
    iE  = [] ;
    iV  = [] ;
    iP  = [] ;
    iI  = [] ;
    
    
    %   Control Volume elasticity parameters
    volumeMax = [];
    isElastic = [];
    
    
    %   Momentum Cell Volume and Mass
    volMCfrom     = []  ;
    volMCto       = []  ;
    volMC         = []  ;
    volMCFracFrom = []  ;
    volMCFracTo   = []  ;
    massMC        = []  ;
    
    %   Conserved quantity sources
    sMass     = [] ;
    sEnergy   = [] ;
    sMomentum = [] ;

    %   Momentum RHS Parameters
    LoD      = [] ;
    flowArea = [] ;
    Ainter   = [] ;
    
    
    %   Various orientation dependent vectors
    upDotN   = []   ;
    downDotN = []   ;
    nX       = []   ;
    nY       = []   ;
    zXup     = []   ;
    zYup     = []   ;
    zXdown   = []   ;
    zYdown   = []   ;
    gizi     = []   ;
    
    
    %   Connectivity summation arrays
    Ccv    = [];
    Cmc    = [];
    Cinter = [];
    
    
    %   Miscellaneous
    time      = 0               ;
    dfdq      = []              ;
    TD        = struct('T',[])  ;
    isDynamic = true()          ;
    
    
    
    
    % =================================================== %
    %                  Method Defintions                  %
    % =================================================== %
    q2Dup.prepare = @(varargin) prepare(varargin{:});
    function [] = prepare(q0,t,varargin)
        if isNotPrepared
            
            %   Load model
            extractModelValues();
            assignSelfValues();
            
            %   Initialize state
            time = t;
            setAll(makeDimensionless(q0));
            updateClosureEnvironment();
            
            isNotPrepared = false;
            
        end
    end
    
    function [] = extractModelValues()
        
        % Conserved quantities
        model    = model.get();
        mass     = model.controlVolume.mass     ;
        energy   = model.controlVolume.energy   ;
        volume   = model.controlVolume.volume   ;
        momentum = model.momentumCell.momentum  ;
        
        %   Pull dimensalizers from model
        massDim     = model.dimensionalizer.mass    ;
        energyDim   = model.dimensionalizer.energy  ;
        volumeDim   = model.dimensionalizer.volume  ;
        momentumDim = model.dimensionalizer.momentum;
        
        %   Momentum cell sense
        from  = model.momentumCell.from ;
        to    = model.momentumCell.to   ;
        
        %   Interface sense
        up    = model.interface.up      ;
        down  = model.interface.down    ;
        
        % Control volume elasticity logical
        volumeMax = model.controlVolume.maximumVolume   ;
        isElastic = volume < volumeMax                  ;
        
        % Volumes of momentum cells
        volMCFracFrom  = model.momentumCell.volumeFractionFrom  ;
        volMCFracTo    = model.momentumCell.volumeFractionTo    ;
        volMCFracTotal = volMCFracFrom + volMCFracTo            ;
        volMCFracFrom  = volMCFracFrom ./ volMCFracTotal        ;
        volMCFracTo    = volMCFracFrom ./ volMCFracTotal        ;
        
        % Interface parameters
        Ainter = model.interface.flowArea   ;
        
        % Sources
        sMass     = model.controlVolume.source.mass     ;
        sEnergy   = model.controlVolume.source.energy   ;
        sMomentum = model.momentumCell.source.momentum  ;
        
        % Momentum
        flowArea = model.momentumCell.flowArea          ;
        LoD      = model.momentumCell.LoD               ;
        
        
        
        
        % ======================================================================= %
        %                            Parameter calculating                        %
        % ======================================================================= %
        
        %   Indices
        nCV = max([from(:);to(:)])  ;
        nMC = length(from)          ;
        iM  = 1:nCV                 ;
        iE  = iM + nCV              ;
        iV  = iE + nCV              ;
        iP  = iV(end) + (1:nMC)          ;
        %
        q2Dup.set('indices.fast',[iM(:);iE(:);iV(:)]);
        q2Dup.set('indices.slow',iP(:));
        
        % (Re-)Normalized of momentum cells
        volMC     = volMCFracFrom .* volume(from) +  volMCFracTo .* volume(to)  ;
        volMCfrom = volMCFracFrom .* volMC                                      ;
        volMCto   = volMCFracTo   .* volMC                                      ;
        
        
        %   Pull orientations
        nX     = model.interface.normalX            ;
        nY     = model.interface.normalY            ;
        zXup   = model.momentumCell.directionX(up)  ;
        zYup   = model.momentumCell.directionY(up)  ;
        zXdown = model.momentumCell.directionX(down);
        zYdown = model.momentumCell.directionY(down);
        
        % Momentum cell-Interface dots
        upDotN   = zXup   .* nX + zYup   .* nY  ;
        downDotN = zXdown .* nX + zYdown .* nY  ;
        
        
        % Gravity
        g     = -9.81;
        gizi  = g * model.momentumCell.directionY;
        
        
        % Summation matrices
        [Ccv,Cmc,Cinter,iI] = ...
            IntrepidTwilight.AdamantWave.toolbox.GetSummationMatrices(...
                [from,to],[up,down],[upDotN,downDotN]);
        
        
        % Initialization for inclusion into the closure environment
        massMC = 0;
        vCV    = 0;
        vI     = 0;
        dfdq   = {zeros(nCV) ; zeros(nCV)  ; zeros(nCV) ; zeros(nMC)};
        
    end
    
    function [] = assignSelfValues()
        epsilon   = q2Dup.get('epsilon')        ;
        isDynamic = q2Dup.get('isDynamicVolume');
        q2Dup.set('isDynamic',...
            [isDynamic;isDynamic;isDynamic.*isElastic;ones(nMC,1)]);
    end
    
    
    
    
    % =================================================== %
    %                  Getters/Setters                    %
    % =================================================== %
    
    %   Set
    function [] = setMass(massStar)
        mass = massStar .* massDim;
    end
    function [] = setEnergy(energyStar)
        energy = energyStar .* energyDim;
    end
    function [] = setVolume(volumeStar)
        volume = volumeStar .* volumeDim;
    end
    function [] = setMomentum(momentumStar)
        momentum = momentumStar .* momentumDim;
    end
    function [] = setAll(q)
        setMass(q(iM));
        setEnergy(q(iE));
        setVolume(q(iV));
        setMomentum(q(iP));
        updateClosureEnvironment();
    end
    function [] = setFast(q)
        setMass(q(iM));
        setEnergy(q(iE));
        setVolume(q(iV));
        updateThermodynamicState();
    end
    function [] = setSlow(q)
        setMomentum(q);
        updateVelocity();
    end
    %   Set
    function q = getMass()
        q = mass ./ massDim ;
    end
    function q = getEnergy()
        q = energy ./ energyDim;
    end
    function q = getVolume()
        q = volume ./ volumeDim;
    end
    function q = getMomentum()
        q = momentum ./ momentumDim;
    end
    function q = getAll()
        q = [...
            getMass()       ;
            getEnergy()     ;
            getVolume()     ;
            getMomentum()   ];
    end
    function q = getFast()
        q = [...
            getMass()   ;
            getEnergy() ;
            getVolume() ];
    end
    function q = getSlow()
        q = getMomentum();
    end
    
    
    
    
    % =================================================== %
    %                  Dimensionalizers                   %
    % =================================================== %
    
    %   Give dimensions
    function q = makeDimensional(q)
        q = [...
            q(iM) .* massDim     ;
            q(iE) .* energyDim   ;
            q(iV) .* volumeDim   ;
            q(iP) .* momentumDim ];
    end
    function q = makeFastDimensional(q)
        q = [...
            q(iM) .* massDim     ;
            q(iE) .* energyDim   ;
            q(iV) .* volumeDim   ];
    end
    function q = makeSlowDimensional(q)
        q = q .* momentumDim;
    end
    %   Take dimensions
    function q = makeDimensionless(q)
        q = [...
            q(iM) ./ massDim     ;
            q(iE) ./ energyDim   ;
            q(iV) ./ volumeDim   ;
            q(iP) ./ momentumDim ];
    end
    function q = makeFastDimensionless(q)
        q = [...
            q(iM) ./ massDim     ;
            q(iE) ./ energyDim   ;
            q(iV) ./ volumeDim   ];
    end
    function q = makeSlowDimensionless(q)
        q = q ./ momentumDim;
    end
    
    
    
    % =================================================== %
    %                     Full RHS                        %
    % =================================================== %
    function f = rhs(q,t)
        
        % Pull conserved values
        mass     = q(iM) .* massDim     ;
        energy   = q(iE) .* energyDim   ;
        volume   = q(iV) .* volumeDim   ;
        momentum = q(iP) .* momentumDim ;
        time     = t                    ;
        
        updateClosureEnvironment();
        
        f = [...
            rhsMass()     ./ massDim   .* isDynamic ;
            rhsEnergy()   ./ energyDim .* isDynamic ;
            rhsVolume()   ./ volumeDim .* isDynamic ;
            rhsMomentum() ./ momentumDim            ];
    end
    
    
    % =================================================== %
    %             Combination Mass-Energy RHS             %
    % =================================================== %
    function f = rhsFast(q,t)
        % Pull conserved values
        mass   = q(iM) .* massDim   ;
        energy = q(iE) .* energyDim ;
        volume = q(iV) .* volumeDim ;
        time   = t                  ;
        
        updateThermodynamicState();
        
        f = [...
            rhsMass()   ./ massDim   .* isDynamic ;
            rhsEnergy() ./ energyDim .* isDynamic ;
            rhsVolume() ./ volumeDim .* isDynamic ];
        
    end
    function f = rhsSlow(q,t)
        
        %	Pull conserved values
        momentum = q .* momentumDim ;
        time     = t                ;
        
        %   Get flux term
        updateVelocity();
        f = rhsMomentum()./momentumDim;
        
    end
    
    
    
    %{
    ===========================================================
                             Mass RHS
    ===========================================================
    %}
    function f = rhsMass(massStar,t)
        
        if (nargin >= 1)
            mass = massStar .* massDim;
        end
        if (nargin >= 2)
            time = t;
        end
        
        % Advection term
        rhod  = (vCV >= 0) .* TD.rho(from)  +...
                (vCV <  0) .* TD.rho(to)    ;
        vzRho =  vCV .* flowArea .* rhod    ;
        
        % Total RHS
        f = Ccv*(vzRho) + sMass(mass,energy,momentum,TD,time) ;
        
    end
    
    
    
    %{
    ===========================================================
                            Energy RHS
    ===========================================================
    %}
    function f = rhsEnergy(energyStar,t)
        
        if (nargin >= 1)
            energy = energyStar .* energyDim;
        end
        if (nargin >= 2)
            time = t;
        end
        
        
        % Advection term
        rhohd  = (vCV >= 0) .* TD.rhoh(from)    +...
                 (vCV <  0) .* TD.rhoh(to)      ;
        vzRhoh =  vCV .* flowArea .* rhohd      ;
        
        % Total RHS
        f  = Ccv*vzRhoh + sEnergy(mass,energy,momentum,TD,time) ;
        
    end
    
    
    
    %{
    ===========================================================
                             Mass RHS
    ===========================================================
    %}
    function f = rhsVolume(volumeStar,~)
        
        if (nargin >= 1)
            volume = volumeStar .* massDim;
        end
        %         if (nargin >= 2)
        %             time = t;
        %         end
        
        % Advection term
        vzVolume = vCV .* flowArea ;
        
        % Total RHS
        f = (Ccv*(vzVolume)) .* isElastic;
        
    end
    
    
    
    %{
    ===========================================================
                          Momentum RHS
    ===========================================================
    %}
    function f = rhsMomentum(momentumStar,t)
        
        if (nargin >= 1)
            momentum = momentumStar .* momentumDim;
        end
        if (nargin >= 2)
            time = t;
        end

        %   Intensive momentum
        rhov = momentum ./ volMC ;
        
        %   Calculate momentum exchange term
        isDownwind = (vI > 0)                                               ;
        isUpwind   = not(isDownwind)                                        ;
        vIrhov     = vI .* (isDownwind .* rhov(up) + isUpwind .* rhov(down));
        
        % Get the friction factor
        args      = FilterList(from,abs(rhov),abs(vCV),1E-4)                            ;
        fFricFrom = frictionFactor(args{:},volMCfrom./flowArea,structFilter(TD,from))   ;
        args      = FilterList(to,abs(rhov),abs(vCV),1E-4)                              ;
        fFricTo   = frictionFactor(args{:},volMCto./flowArea,structFilter(TD,to))       ;
        fFriction = volMCFracFrom .* fFricFrom + volMCFracTo .* fFricTo                 ;
        
        % Pieces
        advect = Cmc*(vIrhov.*Ainter) - Cinter*(TD.P(iI).*Ainter)       ;
        buoy   = gizi .* massMC                                         ;
        fric   = -fFriction .* LoD .* abs(rhov) .* vCV .* flowArea / 2  ;
        
        % Total RHS
        f = advect + buoy + fric + sMomentum(mass,energy,momentum,TD,time) ;
        
    end
    
    
    
    
    
    
    
    %{
    ===========================================================
                          Update functions
    ===========================================================
    %}
    
    function [] = updateClosureEnvironment()
        updateThermodynamicState();
        updateVelocity();
    end
    
    function [] = updateVelocity()

        % Volume-average density and CV surface velocities
        volMCfrom = volMCFracFrom .* volume(from)                           ;
        volMCto   = volMCFracTo   .* volume(to)                             ;
        volMC     = volMCfrom +  volMCto                                    ;
        massMC    = volMCFracFrom .* mass(from) + volMCFracTo .* mass(to)   ;
        vCV       = momentum ./ massMC                                      ;
        
        %   Momentum cell interface velocity
        vI = (momentum(up) .* upDotN + momentum(down) .* downDotN) ./   ...
                         (massMC(up) + massMC(down))                    ;
        
    end
    
    function [] = updateThermodynamicState()
        % Thermodynamic properites
        TD.rho    = mass   ./ volume                ;
        TD.i      = energy ./ mass                  ;
        [~,TD]    = Temperature(TD.rho,TD.i,400+TD.i*0)        ;
        TD.P      = Pressure(TD.rho,TD.T,true,TD)   ;
        TD.rhoh   = energy ./ volume + TD.P         ;
        
        
        
        if (numel(TD.P) <= nCV)
            if any(TD.T > 1.5*TD.T(1))
                g = [];
            end
        end
        
        if any(isnan([TD.P;TD.T]))
            g = [];
        end
        
    end
    
    
    
    %{
    ===========================================================
                          Friction Factor
    ===========================================================
    %}
    
    %   Phase-safe friction factor function
    function f = frictionFactor(rhov,v,roughness,Lchar,TD)
        
        %   Scalar expansion
        if isscalar(Lchar)
            Lchar = ones(size(rhov))*Lchar;
        end
        if isscalar(roughness)
            roughness = ones(size(rhov))*roughness;
        end
        
        %   Allocation
        f = rhov;
        
        %   Phase mask
        isTwoPhi = TD.isTwoPhi  ;
        isOnePhi = not(isTwoPhi);
        
        %   Get kinematic numbers
        mu  = Viscosity(TD.rho,TD.T,true,TD);
        Re  = rhov .* Lchar ./ mu           ;
        eoD = roughness ./ Lchar            ;


        %   One-phase friction factor
        if any(isOnePhi)
            f(isOnePhi) = kinematicFrictionFactor(Re(isOnePhi),eoD(isOnePhi)) ;
        end
        
        %   Two-phase friction factor
        if any(isTwoPhi)
            
            % Friedel correlation
            %   As presented in "Convective Boiling and Condensation" by Collier and
            %   Thome.
            
            %   Mask properties for from
            rhov    = rhov(isTwoPhi)            ;
            v       = v(isTwoPhi)               ;
            Lchar   = Lchar(isTwoPhi)           ;
            eoD     = eoD(isTwoPhi)             ;
            TD      = structFilter(TD,isTwoPhi) ;
            rhoL    = TD.rhoL                   ;
            rhoG    = TD.rhoG                   ;
            x       = TD.x                      ;
            T       = TD.T                      ;
            
            %   Liquid/Gas Reynolds numbers and friction factors
            muL  = Viscosity(rhoL,T,true,TD)            ;
            muG  = Viscosity(rhoG,T,true,TD)            ;
            ReL  = rhov .* Lchar ./ muL                 ;
            ReG  = rhov .* Lchar ./ muG                 ;
            fL0  = kinematicFrictionFactor(ReL,eoD)     ;
            fG0  = kinematicFrictionFactor(ReG,eoD)     ;
             
            %   Correlation for two-phase friction multiplier
            A1     = (1-x).^2 + x.^2 .* ( rhoL.*fG0 ./ (rhoG.*fL0) );
            A2     = x.^0.78 .* (1-x).^0.24;
            A3     = (rhoL./rhoG).^0.91 .* (muG./muL).^0.19 .* (1-muG./muL).^0.7;
            Fr     = v.^2 ./(9.81*Lchar);
            sigma  = SurfaceTension(T);
            We     = (rhov .* v) .* Lchar ./ sigma;
            phi2f0 = A1 + 3.24 * A2.*A3./(Fr.^0.045.*We.^0.035);
            
            %   Two-phase friction factor
            f(isTwoPhi) = phi2f0 .* fL0;
            
        end
        
        f = real(f);
        
    end


    %   Friction factor determined by kinematics alone
    function f = kinematicFrictionFactor(Re,eoD)
        %{

        Implements the smooth laminar equation and turbulent Serghides
        correlation (1984) for the Darcy friction factor as presented in 
        "Explicit Friction Factor Accuracy and Computational Efficiency for 
        Turbulent Flow in Pipes" (2013) by Coole, Winning.

        Log-Linear interpolation is used in the kinematic transition regime.

        %}
        
        %   Initialize
        f = Re;
        
        %   Masks
        isLaminar    = Re <= 2000                   ;
        isTurbulent  = Re >= 4000                   ;
        isTransition = not(isLaminar | isTurbulent) ;
        
        %   Laminar
        if any(isLaminar)
            f(isLaminar) = 64 ./ Re(isLaminar);
        end
        
        %   Turbulent
        if any(isTurbulent)
            ReT = Re(isTurbulent);
            eoDT = eoD(isTurbulent)/3.7;
            a = -2*log(12     ./ ReT + eoDT);
            b = -2*log(2.51*a ./ ReT + eoDT);
            c = -2*log(2.51*b ./ ReT + eoDT);
            f(isTurbulent) = ( a  -  (b-a).^2./(c-2*b+a) ).^(-2);
        end
        
        %   Transition
        if any(isTransition)
            
            %   Laminar side
            fLam = 64./2000;
            
            % Turbulent side
            eoDT = eoD(isTransition)/3.7;
            a     = -2*log10(12     ./ 4000 + eoDT) ;
            b     = -2*log10(2.51*a ./ 4000 + eoDT) ;
            c     = -2*log10(2.51*b ./ 4000 + eoDT) ;
            fTurb = (a-(b-a).^2./(c-2*b+a)).^(-2)   ;
            
            %   Log-linear interpolation
            logRe = log10(Re(isTransition))   ;
            log4t = log10(4000)               ;
            log2t = log10(2000)               ;
            wTurb = (logRe-log2t)/(log4t-log2t);
            wLam  = (log4t-logRe)/(log4t-log2t);
            
            %   Transition factor
            f(isTransition) = wTurb.*fTurb + wLam.*fLam;
            
        end
        
        f = real(f);
        
    end
    
    
    
    
    
    %{
    ===========================================================
                          Jacobians
    ===========================================================
    %}
    
    function dfdqOut = jacobianBlockDiagonal(q)
        
        %   Pull state values
        q           = makeDimensional(q);
        massRef     = q(iM)             ;
        energyRef   = q(iE)             ;
        volumeRef   = q(iV)             ;
        momentumRef = q(iP)             ;
        
        
        %   Perturbed values
        massStep     = massRef     * epsilon;
        energyStep   = energyRef   * epsilon;
        volumeStep   = volumeRef   * epsilon;
        momentumStep = momentumRef * epsilon;
        %
        massTD     = massRef     + massStep     ;
        energyTD   = energyRef   + energyStep   ;
        volumeTD   = volumeRef   + volumeStep   ;
        momentumTD = momentumRef + momentumStep ;
        %
        %   Build large perturbed arrays for fast Thermodynamic update
        mass     = [massRef     ; massTD      ; massRef     ; massRef     ; massRef     ];
        energy   = [energyRef   ; energyRef   ; energyTD    ; energyRef   ; energyRef   ];
        volume   = [volumeRef   ; volumeRef   ; volumeRef   ; volumeTD    ; volumeRef   ];
        momentum = [momentumRef ; momentumRef ; momentumRef ; momentumRef ; momentumTD  ];
        
        
        %   Update and store perturbed values
        TD.T = repmat(TD.T,5,1);
        updateThermodynamicState();
        rhoBlock      = TD.rho      ;
        rhoLBlock     = TD.rhoL     ;
        rhoGBlock     = TD.rhoG     ;
        eBlock        = TD.i        ;
        Tblock        = TD.T        ;
        Pblock        = TD.P        ;
        rhohBlock     = TD.rhoh     ;
        xBlock        = TD.x        ;
        isTwoPhiBlock = TD.isTwoPhi ;
        
        %   Create unmutated TD struct
        TDref = structFilter(TD,1:nCV);
%         TDref.rho      = rhoBlock(1:nCV)        ;
%         TDref.rhoL     = rhoLBlock(1:nCV)       ;
%         TDref.rhoG     = rhoGBlock(1:nCV)       ;
%         TDref.e        = eBlock(1:nCV)          ;
%         TDref.T        = Tblock(1:nCV)          ;
%         TDref.P        = Pblock(1:nCV)          ;
%         TDref.rhoh     = rhohBlock(1:nCV)       ;
%         TDref.x        = xBlock(1:nCV)          ;
%         TDref.isTwoPhi = isTwoPhiBlock(1:nCV)   ;
        
        %   Store compressed perturbed vectors and reset closure variables
        massTD     = [ massRef     ; massTD     ]   ;
        energyTD   = [ energyRef   ; energyTD   ]   ;
        volumeTD   = [ volumeRef   ; volumeTD   ]   ;
        momentumTD = [ momentumRef ; momentumTD ]   ;
        mass       = massRef                        ;
        energy     = energyRef                      ;
        volume     = volumeRef                      ;
        momentum   = momentumRef                    ;
        TD         = TDref                          ;
        
        
        
        
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
            nCVk           = nCV+k;
            mass(k)        = massTD(nCVk)       ;
            TD.rho(k)      = rhoBlock(nCVk)     ;
            TD.rhoL(k)     = rhoLBlock(nCVk)    ;
            TD.rhoG(k)     = rhoGBlock(nCVk)    ;
            TD.e(k)        = eBlock(nCVk)       ;
            TD.T(k)        = Tblock(nCVk)       ;
            TD.P(k)        = Pblock(nCVk)       ;
            TD.rhoh(k)     = rhohBlock(nCVk)    ;
            TD.x(k)        = xBlock(nCVk)       ;
            TD.isTwoPhi(k) = isTwoPhiBlock(nCVk);
            updateVelocity();
            massPlus  = rhsMass();
            
            %   Calculate
            dfdq{1}(:,k) = (massPlus - mass0)/(massStep(k));
            
            %   Reset
            mass(k)    = massRef(k)     ;
            TD.rho(k)  = TDref.rho(k)   ;
            TD.i(k)    = TDref.i(k)     ;
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
            nCVk           = 2*nCV + k          ;
            energy(k)      = energyTD(nCV + k)  ;
            TD.rho(k)      = rhoBlock(nCVk)     ;
            TD.rhoL(k)     = rhoLBlock(nCVk)    ;
            TD.rhoG(k)     = rhoGBlock(nCVk)    ;
            TD.e(k)        = eBlock(nCVk)       ;
            TD.T(k)        = Tblock(nCVk)       ;
            TD.P(k)        = Pblock(nCVk)       ;
            TD.rhoh(k)     = rhohBlock(nCVk)    ;
            TD.x(k)        = xBlock(nCVk)       ;
            TD.isTwoPhi(k) = isTwoPhiBlock(nCVk);
            updateVelocity();
            energyPlus = rhsEnergy();
            
            %   Calculate
            dfdq{2}(:,k) = (energyPlus - energy0)/(energyStep(k));
            
            %   Reset
            energy(k)  = energyRef(k)   ;
            TD.rho(k)  = TDref.rho(k)   ;
            TD.i(k)    = TDref.i(k)     ;
            TD.T(k)    = TDref.T(k)     ;
            TD.P(k)    = TDref.P(k)     ;
            TD.rhoh(k) = TDref.rhoh(k)  ;
            
        end
        
        
        % ========================================================= %
        %                       Volume Block                        %
        % ========================================================= %
        
        %   Get unperturbed value
        updateVelocity();
        volume0 = rhsVolume();
        
        %   Iterate through columns of the block
        for k = K(isElastic)
            
            %   Perturb
            nCVk           = 3*nCV + k          ;
            volume(k)      = volumeTD(nCV+k)    ;
            TD.rho(k)      = rhoBlock(nCVk)     ;
            TD.rhoL(k)     = rhoLBlock(nCVk)    ;
            TD.rhoG(k)     = rhoGBlock(nCVk)    ;
            TD.e(k)        = eBlock(nCVk)       ;
            TD.T(k)        = Tblock(nCVk)       ;
            TD.P(k)        = Pblock(nCVk)       ;
            TD.rhoh(k)     = rhohBlock(nCVk)    ;
            TD.x(k)        = xBlock(nCVk)       ;
            TD.isTwoPhi(k) = isTwoPhiBlock(nCVk);
            updateVelocity();
            volumePlus = rhsVolume();
            
            %   Calculate
            dfdq{3}(:,k) = (volumePlus - volume0)/(volumeStep(k));
            
            %   Reset
            volume(k)  = volumeRef(k)   ;
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
            %             nMCk        = 3*nCV+k;
            momentum(k) = momentumTD(nMC+k) ;
%             TD.rho(k)      = rhoBlock(nCVk)     ;
%             TD.rhoL(k)     = rhoLBlock(nCVk)    ;
%             TD.rhoG(k)     = rhoGBlock(nCVk)    ;
%             TD.e(k)        = eBlock(nCVk)       ;
%             TD.T(k)        = Tblock(nCVk)       ;
%             TD.P(k)        = Pblock(nCVk)       ;
%             TD.rhoh(k)     = rhohBlock(nCVk)    ;
%             TD.x(k)        = xBlock(nCVk)       ;
%             TD.isTwoPhi(k) = isTwoPhiBlock(nCVk);
            updateVelocity();
            momentumPlus = rhsMomentum();
            
            %   Calculate
            dfdq{4}(:,k) = (momentumPlus - momentum0)/(momentumStep(k));
            
            %   Reset
            momentum(k) = momentumRef(k);
%             TD.rho(k)   = TDref.rho(k)  ;
%             TD.i(k)     = TDref.i(k)    ;
%             TD.T(k)     = TDref.T(k)    ;
%             TD.P(k)     = TDref.P(k)    ;
%             TD.rhoh(k)  = TDref.rhoh(k) ;
            
        end
        
        %   Reset all mutations to closure environment
        mass     = massRef      ;
        energy   = energyRef    ;
        volume   = volumeRef    ;
        momentum = momentumRef  ;
        TD       = TDref        ;
        
        
        %   Pass out block diagonal
        dfdqOut = dfdq;
        
    end
    
    
end