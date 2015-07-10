% clc();
% clear();

hem = IntrepidTwilight.AdamantWave.HEM();

r2 = cos(pi/4);
s2 = sqrt(2);


%  12 control volumes
hem.model.momentumCell.from = (1:12)';
hem.model.momentumCell.to   = [2:12,1]';
%
%   Flow direction

hem.model.momentumCell.directionX = [+0;+0;+0;+0;+1;+1;+0;+0;+0;+0;-1;-1];
hem.model.momentumCell.directionY = [-1;-1;-1;-1;+0;+0;+1;+1;+1;+1;+0;+0];
%
%   Interfaces
hem.model.interface.volumeUp       = (1:12)'  ;
hem.model.interface.volumeDown     = [2:12,1]';
hem.model.interface.normalX        = [+0, +0, +0, +r2, +1, +r2, +0, +0, +0, -r2, -1, -r2]';
hem.model.interface.normalY        = [-1, -1, -1, -r2, +0, +r2, +1, +1, +1, +r2, +0, -r2]';
hem.model.interface.surfaceArea    = [+1, +1, +1, +s2, +1, +s2, +1, +1, +1, +s2, +1, +s2]';
hem.model.momentumCell.volumeBack  = ones(12,1)/2;
hem.model.momentumCell.volumeFront = ones(12,1)/2;
hem.model.momentumCell.LoD         = 1;


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
onesCV = ones(problem.miscellaneous.nCV,1) ;
onesMC = ones(problem.miscellaneous.nMC,1) ;

hem.model.controlVolume.volume     = ones(12,1)                             ;
hem.model.controlVolume.mass       = rho0 * hem.model.controlVolume.volume  ;
hem.model.controlVolume.energy     = e0   * hem.model.controlVolume.volume  ;
hem.model.momentumCell.momentum    = v0   * hem.model.controlVolume.volume  ;
hem.model.momentumCell.momentum(1) = 10 * hem.mode.momentumCell.momentum(1) ;

problem.dimensionalizer.rho  = rho0    ;
problem.dimensionalizer.rhoe = rhoe0   ;
problem.dimensionalizer.rhov = rhov0   ;

problem.miscellaneous.iRho  = (1:problem.miscellaneous.nCV)';
problem.miscellaneous.iRhoe = problem.miscellaneous.nCV+ problem.miscellaneous.iRho ;
problem.miscellaneous.iRhov = 2*problem.miscellaneous.nCV + (1:problem.miscellaneous.nMC)';


% Thermodynamic properites
problem.initialState.e = problem.initialState.rhoe0 ./ problem.initialState.rho0     ;
problem.initialState.T = Temperature(problem.initialState.rho0,problem.initialState.e);
problem.initialState.P = Pressure(problem.initialState.rho0,problem.initialState.T)   ;

problem.miscellaneous.epsilon = 1E-8;



%   Method choices
problem.semidiscretization.name         = 'Quasi2DUpwind'   ;
problem.timeStepper.name                = 'implicitEuler'   ;
problem.timeStepper.stepSize            = 0.2               ;
problem.solver.name                     = 'JFNK'            ;
% problem.solver.preconditioner.type      = 'none'            ;
problem.solver.preconditioner.type      = 'block-jacobi'    ;
% problem.solver.preconditioner.type      = 'full-stagnant'   ;
problem.solver.preconditioner.blockSize = [problem.miscellaneous.nCV;problem.miscellaneous.nCV;problem.miscellaneous.nMC];


qHi = [996.8/rho0*onesCV;Inf*onesCV;Inf*onesMC];
qLo = [1E-5*onesCV;-Inf*onesCV;-Inf*onesMC];
problem.solver.guard.value = @(q)    guardValue(q,qLo,qHi,problem);
problem.solver.guard.step  = @(q,dq) guardStep(q,dq,qLo,qHi,problem);

% tic;
% simulation.run([0,1],0.01,0.1);
% toc;



% sd     = IntrepidTwilight.AdamantWave.Quasi2DUpwind(problem);
% dfFull = IntrepidTwilight.ConvenientMeans.numericalJacobianFull(sd.rhs,problem.initialState.q0);
% dfFull = eye(size(dfFull)) - 0.1*dfFull;
% 
% Lambda = eig(dfFull,'vector');
% Show(real(Lambda))
