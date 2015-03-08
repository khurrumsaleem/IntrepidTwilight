clc();
% clear();

problem = IntrepidTwilight.new('problem');

%   Upfront problem choices
problem.semidiscretization.name = 'Quasi2DUpwind' ;
problem.timeStepper.name        = 'implicitEuler' ;
problem.timeStepper.stepSize    = 0.2             ;
problem.solver.name             = 'JFNK'          ;
problem.solver.preconditioner   = 'BlockJacobi'   ;

r2 = cos(pi/4);
s2 = sqrt(2);

% %   19 control volumes
% s.from     = [1:8, 7,10:11,12,9 ,13:18,19]';
% s.to       = [2:9,10,11:12,18,13,14:19,20]';
% s.zx       = [ 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, -1, -1 ]';
% s.zy       = [-1,-1,-1,-1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1,  0,  0,  0,  0 ]';
% %
% s.up       = [1:7               , 6, 7,  8, 9:11,12,12,13:19]';
% s.down     = [2:8               , 9, 9, 13,10:12,18,19,14:20]';
% s.nx       = [ 0, 0, 0, r2,1,1,1,r2,-r2,r2,0,0,0,r2,-r2,0,0,0,-r2,-1,-1,-1]';
% s.ny       = [-1,-1,-1,-r2,0,0,0,r2, r2,r2,1,1,1,r2, r2,1,1,1, r2, 0, 0, 0]';
% %
% s.volFrom  = [0.5,0.5,0.5,0.5,0.5,0.5,3/7,0.5,1/3,0.5,0.5,2/3,0.5,0.5,0.5,0.5,0.5,0.5,3/7,0.5]';
% s.Ainter   = [1,1,1,s2,1,0.5,1,r2,r2,s2,1,1,1,r2,r2,1,1,1,r2,1,0.5,1]';

%  6 control volumes
problem.geometry.from        = [ 1,2,3,4, 5,6, 3]';
problem.geometry.to          = [ 2,3,4,5, 6,1, 6]';
problem.geometry.zx          = [ 0,1,1,0,-1,-1, 0 ]';
problem.geometry.zy          = [-1,0,0,1, 0, 0, 1 ]';
problem.geometry.up          = [  1,  2,   3,   4,   5,   6,   2,   3,   7,   7]';
problem.geometry.down        = [  2,  3,   4,   5,   6,   1,   7,   7,   5,   6]';
problem.geometry.nx          = [+r2,  1, +r2, -r2,  -1, -r2, +r2, -r2, +r2, -r2]';
problem.geometry.ny          = [-r2,  0, +r2, +r2,   0, -r2, +r2, +r2, +r2, +r2]';
problem.geometry.Ainter      = [ s2,1/2, +s2,  s2, 1/2,  s2,  r2,  r2,  r2,  r2]';
problem.geometry.volumeBack  = [1/2,1/2, 3/8, 1/2, 1/2, 3/8, 1/4]';
problem.geometry.volumeFront = [1/2,3/8, 1/2, 1/2, 3/8, 1/2, 1/4]';
problem.geometry.LoD         = 1;


problem.miscellaneous.nCV      = max([problem.geometry.from;problem.geometry.to]);
problem.miscellaneous.nMC      = length(problem.geometry.from);
problem.miscellaneous.nInter   = length(problem.geometry.nx);
problem.miscellaneous.nEq      = 2*problem.miscellaneous.nCV + problem.miscellaneous.nMC;


problem.miscellaneous.friction = 0.01;

problem.miscellaneous.sRho  = 0;
problem.miscellaneous.sRhov = 0;
problem.miscellaneous.sRhoe = 10E6*[-1;0;0;1;0;0];

rho0 = 9.965569351080000e+02;
e0   = 1.125536123942350e+05;
v0   = 0.01;
problem.initialState.rho0     = rho0*ones(problem.miscellaneous.nCV,1)     ;
problem.initialState.rhoe0    = rho0*e0*ones(problem.miscellaneous.nCV,1)  ;
problem.initialState.rhov0    = rho0*v0*ones(problem.miscellaneous.nMC,1)  ;
problem.initialState.rhov0(1) = 10*problem.initialState.rhov0(1)          ;

problem.initialState.q0 = [problem.initialState.rho0;problem.initialState.rhoe0;problem.initialState.rhov0];

problem.miscellaneous.iRho  = (1:problem.miscellaneous.nCV)';
problem.miscellaneous.iRhoe = problem.miscellaneous.nCV+ problem.miscellaneous.iRho ;
problem.miscellaneous.iRhov = (2*problem.miscellaneous.nCV + 1:problem.miscellaneous.nMC)';


% Thermodynamic properites
problem.initialState.e = problem.initialState.rhoe0 ./ problem.initialState.rho0     ;
problem.initialState.T = Temperature(problem.initialState.rho0,problem.initialState.e);
problem.initialState.P = Pressure(problem.initialState.rho0,problem.initialState.T)   ;

problem.miscellaneous.epsilon = 1E-4;

Simulation = IntrepidTwilight.new('simulation',problem);


%{
s.epsilon = 5E-4;

upwind = IntrepidTwilight.AdamantWave.Semidiscretizations.Quasi2DUpwind(s);
df     = upwind.blockDiagonalJacobian(s.q0);

%   Jacobian sensistivity to epsilon 1E-4 or 1E-5 look good.

epsilons = 10.^((-3:-1:-10)')        ;
df       = zeros(length(epsilons),1) ;
for k = 1:length(epsilons)
    s.epsilon = epsilons(k);
    upwind    = IntrepidTwilight.AdamantWave.Semidiscretizations.Quasi2DUpwind(s);
    df(k)     = norm(upwind.blockDiagonalJacobian(s.q0));
end
loglog(epsilons,df);
%}


% r    = @(qND) (qND - q0./qstar) - s.dt*SemidiscreteUpwind(qND.*qstar,s)./qstar;
% r    = @(q) (q - q0) - s.dt*SemidiscreteUpwind(q,s);
% drdq = CalculateJacobian(r,q0,1E-8*q0);


%{
DensityGuard = @(rhoe) rhoe > [996.8/rho0*ones(s.nCV,1);Inf*ones(s.nCV,1);Inf*ones(s.nMC,1)];


loop = 0;
while true
    fprintf('Time = %d:\n',(loop+1)*s.dt);
    qND = JFNKHouseholder(q0./qstar,r,1E-8,DensityGuard);
    
    %   Update guess temperature
    q0    = qND .* qstar        ;
    rho   = q0(s.rhoMask)       ;
    rhoe  = q0(s.rhoeMask)      ;
    s.T   = Temperature(rho,rhoe./rho,s.T);
    
    %   Rebuild residual
    r  = @(qND) (qND - q0./qstar) - s.dt*SemidiscreteUpwind(qND.*qstar,s)./qstar;
    
    loop = loop + 1;
end

%}
