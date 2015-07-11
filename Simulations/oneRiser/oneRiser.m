% clc();
% clear();

hem = IntrepidTwilight.AdamantWave.HEM();

r2 = cos(pi/4);
s2 = sqrt(2);


%  12 control volumes
hem.modelValue('momentumCell','from', (1:12)');
hem.modelValue('momentumCell','to'  , [2:12,1]');
%
%   Flow direction

hem.modelValue('momentumCell','directionX',[+0;+0;+0;+0;+1;+1;+0;+0;+0;+0;-1;-1]);
hem.modelValue('momentumCell','directionY',[-1;-1;-1;-1;+0;+0;+1;+1;+1;+1;+0;+0]);
%
%   Interfaces
hem.modelValue('interface','up', (1:12)'  );
hem.modelValue('interface','down', [2:12,1]');
hem.modelValue('interface','normalX', [+0, +0, +0, +r2, +1, +r2, +0, +0, +0, -r2, -1, -r2]');
hem.modelValue('interface','normalY', [-1, -1, -1, -r2, +0, +r2, +1, +1, +1, +r2, +0, -r2]');
hem.modelValue('interface','flowArea', [+1, +1, +1, +s2, +1, +s2, +1, +1, +1, +s2, +1, +s2]');
hem.modelValue('momentumCell','volumeFrom', ones(12,1)/2);
hem.modelValue('momentumCell','volumeTo', ones(12,1)/2);
hem.modelValue('momentumCell','LoD', 1);
hem.modelValue('momentumCell','source','friction',0.1);

%   Sources
hem.modelValue('controlVolume','source','mass'   , @(varargin) 0);
hem.modelValue('controlVolume','source','energy' , @(varargin) 0);
hem.modelValue('momentumCell','source','momentum', @(varargin) 0);

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
v0   = 1E-5;
rhoe0 = rho0 * e0;
rhov0 = rho0 * v0;

volume      = ones(12,1)        ;
mass        = rho0  * ones(12,1);
energy      = rhoe0 * ones(12,1);
momentum    = rhov0 * ones(12,1);
momentum(1) = 10*momentum(1)  	;
hem.modelValue('controlVolume','volume' , ones(12,1));
hem.modelValue('controlVolume','mass'   , mass      );
hem.modelValue('controlVolume','energy' , energy    );
hem.modelValue('momentumCell','momentum', momentum  );

% hem.modelValue('dimensionalizer','mass')     = rho0 ;
% hem.modelValue('dimensionalizer','energy')   = rhoe0;
% hem.modelValue('dimensionalizer','momentum') = rhov0;

hem.build();

hem.evolverValue('time','span'          , [0,2]                 )   ;
hem.evolverValue('time','step','maximum', 1                     )   ;
hem.evolverValue('time','step','minimum', 1E-7                  )   ;
hem.evolverValue('time','step','goal'   , 0.1                   )   ;
hem.evolverValue('initialCondition'     , {momentum;[mass;energy]})   ;
hem.evolverValue('saveRate'             , 0.1                   )   ;

hem.evolve();

%   Method choices
% problem.semidiscretization.name         = 'Quasi2DUpwind'   ;
% problem.timeStepper.name                = 'implicitEuler'   ;
% problem.timeStepper.stepSize            = 0.2               ;
% problem.solver.name                     = 'JFNK'            ;
% problem.solver.preconditioner.type      = 'none'            ;
% problem.solver.preconditioner.type      = 'block-jacobi'    ;
% problem.solver.preconditioner.type      = 'full-stagnant'   ;
% problem.solver.preconditioner.blockSize = [problem.miscellaneous.nCV;problem.miscellaneous.nCV;problem.miscellaneous.nMC];


% qHi = [996.8/rho0*onesCV;Inf*onesCV;Inf*onesMC];
% qLo = [1E-5*onesCV;-Inf*onesCV;-Inf*onesMC];
% problem.solver.guard.value = @(q)    guardValue(q,qLo,qHi,problem);
% problem.solver.guard.step  = @(q,dq) guardStep(q,dq,qLo,qHi,problem);

% tic;
% simulation.run([0,1],0.01,0.1);
% toc;



% sd     = IntrepidTwilight.AdamantWave.Quasi2DUpwind(problem);
% dfFull = IntrepidTwilight.ConvenientMeans.numericalJacobianFull(sd.rhs,problem.initialState.q0);
% dfFull = eye(size(dfFull)) - 0.1*dfFull;
% 
% Lambda = eig(dfFull,'vector');
% Show(real(Lambda))
