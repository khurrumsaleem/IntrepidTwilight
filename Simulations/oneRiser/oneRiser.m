clc();
clear();

problem = IntrepidTwilight.new('problem');

r2 = cos(pi/4);
s2 = sqrt(2);


%  12 control volumes
problem.geometry.from        = (1:12)';
problem.geometry.to          = [2:12,1]';
%
%   Flow direction
problem.geometry.zx          = [+0;+0;+0;+0;+1;+1;+0;+0;+0;+0;-1;-1];
problem.geometry.zy          = [-1;-1;-1;-1;+0;+0;+1;+1;+1;+1;+0;+0];
%
%   Interfaces
problem.geometry.up          = (1:12)'  ;
problem.geometry.down        = [2:12,1]';
problem.geometry.nx          = [+0, +0, +0, +r2, +1, +r2, +0, +0, +0, -r2, -1, -r2]';
problem.geometry.ny          = [-1, -1, -1, -r2, +0, +r2, +1, +1, +1, +r2, +0, -r2]';
problem.geometry.Ainter      = [+1, +1, +1, +s2, +1, +s2, +1, +1, +1, +s2, +1, +s2]';
problem.geometry.volumeBack  = ones(12,1)/2;
problem.geometry.volumeFront = ones(12,1)/2;
problem.geometry.LoD         = 1;


problem.miscellaneous.nCV      = max([problem.geometry.from;problem.geometry.to]);
problem.miscellaneous.nMC      = length(problem.geometry.from);
problem.miscellaneous.nInter   = length(problem.geometry.nx);
problem.miscellaneous.nEq      = 2*problem.miscellaneous.nCV + problem.miscellaneous.nMC;


problem.miscellaneous.friction = 0.01;

problem.miscellaneous.sRho  = @(rho,rhoe,rhov,TD,t) 0;
problem.miscellaneous.sRhov = @(rho,rhoe,rhov,TD,t) 0;

energySource = (1:12)'*0;
energySource(8) = +5E7;
problem.miscellaneous.sRhoe = @(rho,rhoe,rhov,TD,t) ((t > 0.01)&&(t<=0.02))*(energySource*(t-0.01)/(0.02-0.01)) + (t>0.02)*energySource;

rho0 = 9.965569351080000e+02;
e0   = 1.125536123942350e+05;
v0   = 0.01;
rhoe0 = rho0 * e0;
rhov0 = rho0 * v0;
onesCV = ones(problem.miscellaneous.nCV,1) ;
onesMC = ones(problem.miscellaneous.nMC,1) ;
problem.initialState.rho0     = rho0*onesCV    ;
problem.initialState.rhoe0    = rho0*e0*onesCV  ;
problem.initialState.rhov0    = rho0*v0*onesMC  ;
problem.initialState.rhov0(1) = 10*problem.initialState.rhov0(1)          ;

problem.initialState.q0 = [...
    problem.initialState.rho0/rho0   ;...
    problem.initialState.rhoe0/rhoe0  ;...
    problem.initialState.rhov0/rhov0  ];
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

problem.miscellaneous.epsilon = 1E-4;



%   Method choices
problem.semidiscretization.name         = 'Quasi2DUpwind'   ;
problem.timeStepper.name                = 'implicitEuler'   ;
problem.timeStepper.stepSize            = 0.2               ;
problem.solver.name                     = 'JFNK'            ;
problem.solver.preconditioner.type      = 'none'            ;
problem.solver.preconditioner.blockSize = [problem.miscellaneous.nCV;problem.miscellaneous.nCV;problem.miscellaneous.nMC];

qHi = [996.8/rho0*onesCV;Inf*onesCV;Inf*onesMC];
qLo = [1E-5*onesCV;-Inf*onesCV;-Inf*onesMC];
problem.solver.guard.value = @(q)    guardValue(q,qLo,qHi);
problem.solver.guard.step  = @(q,dq) guardStep(q,dq,qLo,qHi);

simulation = IntrepidTwilight.new('simulation',problem);

simulation.run([0,1],0.05,0.05);


