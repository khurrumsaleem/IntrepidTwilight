% clc();
% clear();

hem = IntrepidTwilight.AdamantWave.HEM();

r2 = cos(pi/4);
s2 = sqrt(2);


%  12 control volumes
hem.set.model('momentumCell.from', (1:12)');
hem.set.model('momentumCell.to'  , [2:12,1]');
%
%   Flow direction

hem.set.model('momentumCell.directionX',[+0;+0;+0;+0;+1;+1;+0;+0;+0;+0;-1;-1]);
hem.set.model('momentumCell.directionY',[-1;-1;-1;-1;+0;+0;+1;+1;+1;+1;+0;+0]);
hem.set.model('momentumCell.volumeFrom', ones(12,1)/2);
hem.set.model('momentumCell.volumeTo', ones(12,1)/2);
hem.set.model('momentumCell.LoD', 1);
hem.set.model('momentumCell.source.friction',0.1);
%
%   Interfaces
hem.set.model('interface.up', (1:12)'  );
hem.set.model('interface.down', [2:12,1]');
hem.set.model('interface.normalX', [+0, +0, +0, +r2, +1, +r2, +0, +0, +0, -r2, -1, -r2]');
hem.set.model('interface.normalY', [-1, -1, -1, -r2, +0, +r2, +1, +1, +1, +r2, +0, -r2]');
hem.set.model('interface.flowArea', [+1, +1, +1, +s2, +1, +s2, +1, +1, +1, +s2, +1, +s2]');

%   Sources
hem.set.model('controlVolume.source.mass'   , @(varargin) 0);
hem.set.model('controlVolume.source.energy' , @(varargin) 0);
hem.set.model('momentumCell.source.momentum', @(varargin) 0);

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
hem.set.model('controlVolume.volume' , volume    );
hem.set.model('controlVolume.mass'   , mass      );
hem.set.model('controlVolume.energy' , energy    );
hem.set.model('momentumCell.momentum', momentum  );

% hem.set.model('dimensionalizer','mass')     = rho0 ;
% hem.set.model('dimensionalizer','energy')   = rhoe0;
% hem.set.model('dimensionalizer','momentum') = rhov0;

hem.build();

hem.set.evolver('time.span'        , [0,2]                   )   ;
hem.set.evolver('time.step.maximum', 1                       )   ;
hem.set.evolver('time.step.minimum', 1E-7                    )   ;
hem.set.evolver('time.step.goal'   , 0.1                     )   ;
hem.set.evolver('initialCondition' , {[mass;energy];momentum})   ;
hem.set.evolver('saveRate'         , 0.1                     )   ;

hem.evolve();


