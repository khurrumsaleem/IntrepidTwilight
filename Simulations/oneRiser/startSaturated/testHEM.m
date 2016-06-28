function hem = testHEM()
    %   Numbers for normals
    r2 = cos(pi/4);
    s2 = sqrt(2);
    
    %   Simple geometry stuff
    hcorn   = 0.1       ;
    hhorz   = 0.2       ;
    hvert   = 0.3       ;
    Af      = 0.1^2     ;
    volCorn = hcorn * Af;
    volHorz = hhorz * Af;
    volVert = hvert * Af;
    
    %   Vertical / horizontal maps
    vertMap = [-1;-1;-1;-1;-1; 0; 0;+1;+1;+1;+1;+1; 0; 0];
    horzMap = [ 0; 0; 0; 0; 0;+1;+1; 0; 0; 0; 0; 0;-1;-1];
    
    %   Path lengths / hydraulic diameters
    L  = [   hcorn/2 + hvert/2   ;
        hvert*ones(3,1)     ;
        hcorn/2 + hvert/2   ];
    L  = [L ; (hcorn + hhorz)*[1;1]/2 ; L ; (hcorn + hhorz)*[1;1]/2 ];
    Dh = 0.1 * ones(14,1);
    
    %   Define initial state
    T0   = 373.15                     ; % [K]
    x0   = 0.01;
    [P0,rhoL,rhoG] = SaturationStateGivenTemperature(T0);
    rho0 = 1./MixtureProperty(x0,1./rhoL,1./rhoG);
    iL   = InternalEnergyOne(rhoL,T0);
    iG   = InternalEnergyOne(rhoG,T0);
    i0   = MixtureProperty(x0,iL,iG) ;
    v0   = 0.1                       ; % [m/s]
    
    %   Get densities and energies assuming incompressible, hydrostatic pressures
    P     = P0*ones(14,1)   ;
    rho   = P*0 + rho0      ;
    i     = P*0 + i0        ;
    
    
    %   Define extensive properties
    volume       = [              volCorn           ;
        volVert*ones(4,1)      ;
        volCorn ; volHorz ; volCorn ;
        volVert*ones(4,1)      ;
        volCorn ; volHorz      ];
    volumeMax    = volume                                                   ;
    mass         = rho .* volume                                            ;
    energy       = i   .* mass                                              ;
    momentum     = mass .* v0                                               ;
    
    
    %   Instantiate the simulation
    hem = IntrepidTwilight.AdamantWave.HEM();
    
    %  12 control volumes
    hem.set('model','momentumCell.from', (1:14).');
    hem.set('model','momentumCell.to'  , [2:14,1].');
    %
    %   momentum cell geometry
    hem.set('model','momentumCell.directionX'         , horzMap         );
    hem.set('model','momentumCell.directionY'         , vertMap         );
    hem.set('model','momentumCell.volumeFractionFrom' , ones(14,1)/2    );
    hem.set('model','momentumCell.volumeFractionTo'   , ones(14,1)/2    );
    hem.set('model','momentumCell.flowArea'           , Af*ones(14,1)   );
    hem.set('model','momentumCell.LoD'                , L./Dh           );
    hem.set('model','momentumCell.source.friction'    , 0.1             );
    %
    %   Interfaces
    hem.set('model','interface.up'  , (1:14).'  );
    hem.set('model','interface.down', [2:14,1].');
    hem.set('model','interface.normalX',       [  0 ;  0 ;  0 ;  0 ; +r2 ; +1 ; +r2 ;  0 ;  0 ;  0 ;  0 ; -r2 ; -1 ; -r2]);
    hem.set('model','interface.normalY',       [ -1 ; -1 ; -1 ; -1 ; -r2 ;  0 ; +r2 ; +1 ; +1 ; +1 ; +1 ; +r2 ;  0 ; -r2]);
    hem.set('model','interface.flowArea', Af * [ +1 ; +1 ; +1 ; +1 ; +s2 ; +1 ; +s2 ; +1 ; +1 ; +1 ; +1 ; +s2 ; +1 ; +s2]);
    
    %   Sources
    hem.set('model','controlVolume.source.mass'   , @(varargin) 0);
    tramp  = 50;
    shape  = @(t) (t<=tramp)*(t>=10)*(t-10)/(tramp-10) + (t>tramp);
    source = @(t) 5E3 * shape(t);
    hem.set('model','controlVolume.source.energy' , @(~,~,~,~,t) ...
        [0; -source(t) ; zeros(6,1) ; source(t) ; zeros(5,1)]);
    hem.set('model','momentumCell.source.momentum', @(varargin)0);
%         @(~,~,momentum,TD,~)...
%         [0;sign(momentum(2))*Af*(TD.P(2) - P0) ; zeros(12,1)] ...
%         );
    
    %   State
    hem.set('model','controlVolume.mass'          , mass      );
    hem.set('model','controlVolume.energy'        , energy    );
    hem.set('model','controlVolume.volume'        , volume    );
    hem.set('model','controlVolume.maximumVolume' , volumeMax );
    hem.set('model','momentumCell.momentum'       , momentum  );
    
    %   Non-dimensionalizers
    hem.set('model','dimensionalizer.mass'     , max(mass)      );
    hem.set('model','dimensionalizer.energy'   , max(energy)    );
    hem.set('model','dimensionalizer.volume'   , max(volume)    );
    hem.set('model','dimensionalizer.momentum' , max(momentum)  );
    
    %   Set guard
    hem.set('residual','guard.step', @(q,dq) guardStep(q,dq));
    
    %   Preconditioning
    hem.set('preconditioner','kind','block-jacobi');
    
    %   Solver
    hem.set('solver','tolerance.residual', 1E-6);
    hem.set('solver','maximumIterations' , 10);
   
    
    %   Set evolution parameters
    hem.set('evolver','initialCondition'  , [mass;energy;volume;momentum]   );
    hem.set('evolver','tolerance.relative',0.10);
    hem.set('evolver','tolerance.absolute',100 );
end


function dq = guardStep(q,dq)
    while any( (q(1:28) - dq(1:28)) < 0 )
        dq = 0.5 * dq;
    end
end




