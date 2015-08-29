% clc();
% clear();

hem = IntrepidTwilight.AdamantWave.HEM();

r2 = cos(pi/4);
s2 = sqrt(2);


%  12 control volumes
hem.model.set('momentumCell.from', (1:12)');
hem.model.set('momentumCell.to'  , [2:12,1]');
%
%   Flow direction

hem.model.set('momentumCell.directionX',[+0;+0;+0;+0;+1;+1;+0;+0;+0;+0;-1;-1]);
hem.model.set('momentumCell.directionY',[-1;-1;-1;-1;+0;+0;+1;+1;+1;+1;+0;+0]);
hem.model.set('momentumCell.volumeFrom', ones(12,1)/2);
hem.model.set('momentumCell.volumeTo', ones(12,1)/2);
hem.model.set('momentumCell.LoD', 1);
hem.model.set('momentumCell.source.friction',0.1);
%
%   Interfaces
hem.model.set('interface.up', (1:12)'  );
hem.model.set('interface.down', [2:12,1]');
hem.model.set('interface.normalX', [+0, +0, +0, +r2, +1, +r2, +0, +0, +0, -r2, -1, -r2]');
hem.model.set('interface.normalY', [-1, -1, -1, -r2, +0, +r2, +1, +1, +1, +r2, +0, -r2]');
hem.model.set('interface.flowArea', [+1, +1, +1, +s2, +1, +s2, +1, +1, +1, +s2, +1, +s2]');

%   Sources
hem.model.set('controlVolume.source.mass'   , @(varargin) 0);
hem.model.set('controlVolume.source.energy' , @(varargin) 0);
hem.model.set('momentumCell.source.momentum', @(varargin) 0);

% problem.miscellaneous.nCV      = max([problem.geometry.from;problem.geometry.to]);
% problem.miscellaneous.nMC      = length(problem.geometry.from);
% problem.miscellaneous.nInter   = length(problem.geometry.nx);
% problem.miscellaneous.nEq      = 2*problem.miscellaneous.nCV + problem.miscellaneous.nMC;

% problem.miscellaneous.sRho  = @(rho,rhoe,rhov,TD,t) 0;
% problem.miscellaneous.sRhov = @(rho,rhoe,rhov,TD,t) 0;

% energySource = (1:12)'*0;
% energySource(8) = +5E7;
% problem.miscellaneous.sRhoe = @(rho,rhoe,rhov,TD,t) 0;%((t > 0.01)&&(t<=0.02))*(energySource*(t-0.01)/(0.02-0.01)) + (t>0.02)*energySource;

rho0 = 9.965569351080000e+02;
e0   = 1.125536123942350e+05;
v0   = 0.1;
rhoe0 = rho0 * e0;
rhov0 = rho0 * v0;

volume      = 0.1  .* ones(12,1);
mass        = rho0 .* volume    ;
energy      = e0   .* mass      ;
momentum    = mass .* v0        ;
momentum(1) = 10*momentum(1)  	;
hem.model.set('controlVolume.volume' , volume    );
hem.model.set('controlVolume.mass'   , mass      );
hem.model.set('controlVolume.energy' , energy    );
hem.model.set('momentumCell.momentum', momentum  );

% hem.model.set('dimensionalizer','mass')     = rho0 ;
% hem.model.set('dimensionalizer','energy')   = rhoe0;
% hem.model.set('dimensionalizer','momentum') = rhov0;

hem.build();

hem.evolver.set('time.span'        , [0,2]                   )   ;
hem.evolver.set('time.step.maximum', 1                       )   ;
hem.evolver.set('time.step.minimum', 1E-7                    )   ;
hem.evolver.set('time.step.goal'   , 0.1                     )   ;
hem.evolver.set('initialCondition' , {[mass;energy];momentum})   ;
hem.evolver.set('saveRate'         , 0.1                     )   ;

hem.evolve();


