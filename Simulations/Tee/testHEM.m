function hem = testHEM(T1,P1,P7,P10)
    %   Numbers for normals
    r2 = cos(pi/4);
    
    %   Simple geometry stuff
    h   = 0.1   ;
    Af  = h^2   ;
    vol = h * Af;
    
    %   Vertical / horizontal maps
    horzMap = [  ones(6,1) ; zeros(3,1) ];
    vertMap = [ zeros(6,1) ;  ones(3,1) ];
    
    %   Path lengths / hydraulic diameters
    L                  = h*ones(9,1)                                    ;
    Dh                 = h*ones(9,1)                                    ;
    volumeFractionFrom = [ 1 ; 1 ;  1  ; 3/4 ; 1 ; 1 ; 1/2 ; 1 ; 1 ]/2  ;
    volumeFractionTo   = [ 1 ; 1 ; 3/4 ;  1  ; 1 ; 1 ;  1  ; 1 ; 1 ]/2  ;
    volume             = vol * ones(10,1)                               ;
    
    %   State
    v0   = 0.1                      ; % [m/s]
    P    = linspace(P1,P7,7).'      ;
    Ptmp = linspace(P(4),P10,4).'   ;
    P    = [P;Ptmp(2:end)]          ;
    rho  = Density(P,P*0 + T1,P*0)  ;
    i    = InternalEnergy(rho,T1)   ;

    %   Define extensive properties
    volumeMax    = volume           ;
    mass         = rho  .* volume   ;
    energy       = i    .* mass     ;
    momentum     = mass .* v0       ;
    
    
    %   Instantiate the simulation
    hem = IntrepidTwilight.AdamantWave.HEM();
    
    %  12 control volumes
    hem.set('model','momentumCell.from', [1:6,4,8,9].');
    hem.set('model','momentumCell.to'  , [2:7,8,9,10].');
    %
    %   momentum cell geometry
    hem.set('model','momentumCell.directionX'         , horzMap             );
    hem.set('model','momentumCell.directionY'         , vertMap             );
    hem.set('model','momentumCell.volumeFractionFrom' , volumeFractionFrom  );
    hem.set('model','momentumCell.volumeFractionTo'   , volumeFractionTo    );
    hem.set('model','momentumCell.flowArea'           , Af*ones(9,1)        );
    hem.set('model','momentumCell.LoD'                , L./Dh               );

    %
    %   Interfaces
    hem.set('model','interface.up'  , [1:5,3,4,7,8].');
    hem.set('model','interface.down', [2:6,7,7,8,9].');
    hem.set('model','interface.normalX', [ +ones(5,1) ; +r2 ; -r2 ;  0 ;  0 ]);
    hem.set('model','interface.normalY', [ +ones(5,1) ; +r2 ; +r2 ; +1 ; +1 ]);
    hem.set('model','interface.flowArea', Af * [ 1 ; 1 ; 1/2 ; 1 ; 1 ; r2 ; r2 ; 1 ; 1 ]);
    
    %   Sources
    hem.set('model','controlVolume.source.mass'   , @(varargin) 0);
    hem.set('model','controlVolume.source.energy' , @(varargin) 0);
    hem.set('model','momentumCell.source.momentum', @(varargin) 0);
    
    %   State
    hem.set('model','controlVolume.mass'          , mass     );
    hem.set('model','controlVolume.energy'        , energy   );
    hem.set('model','controlVolume.volume'        , volume    );
    hem.set('model','controlVolume.maximumVolume' , volumeMax );
    hem.set('model','momentumCell.momentum'       , momentum );
    
    %   Non-dimensionalizers
    hem.set('model','dimensionalizer.mass'     , mass(:)          );
    hem.set('model','dimensionalizer.energy'   , energy(:)        );
    hem.set('model','dimensionalizer.volume'   , volume(:)        );
    hem.set('model','dimensionalizer.momentum' , 1          );
    
    %
    hem.set('semidiscretization','isDynamicVolume',[false;true(5,1);false;true;true;false]);
    
    %   Set guard
    hem.set('residual','guard.step', @(q,dq) guardStep(q,dq));
    
    %   Preconditioning
    hem.set('preconditioner','kind','block-jacobi');
    
    %   Solver
    hem.set('solver','isNormNotDone',@(x,r) ...
        (norm(r(1:14*3),2)      > 1E-6 ) || ... % Control volume relative measure
        norm(r((14*3+1):end),1) > 1E-4 );
    hem.set('solver','maximumIterations' , 10);
   
    
    %   Set evolution parameters
    hem.set('evolver','initialCondition'  , [mass;energy;volume;momentum]   );
    hem.set('evolver','tolerance.relative',0.10);
    hem.set('evolver','tolerance.absolute',100 );
end

function dq = guardStep(q,dq)
    while any( (q(1:10*3) - dq(1:10*3)) < 0 )
        dq = 0.5 * dq;
    end
end




