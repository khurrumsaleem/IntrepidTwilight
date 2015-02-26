clc;
clear;

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
s.from     = [ 1,2,3,4, 5,6, 3]';
s.to       = [ 2,3,4,5, 6,1, 6]';
s.zx       = [ 0,1,1,0,-1,-1, 0 ]';
s.zy       = [-1,0,0,1, 0, 0, 1 ]';
s.up       = [  1,  2,   3,   4,   5,   6,   2,   3,   7,   7]';
s.down     = [  2,  3,   4,   5,   6,   1,   7,   7,   5,   6]';
s.nx       = [+r2,  1, +r2, -r2,  -1, -r2, +r2, -r2, +r2, -r2]';
s.ny       = [-r2,  0, +r2, +r2,   0, -r2, +r2, +r2, +r2, +r2]';
s.Ainter   = [ s2,1/2, +s2,  s2, 1/2,  s2,  r2,  r2,  r2,  r2]';
s.back     = [1/2,1/2, 3/8, 1/2, 1/2, 3/8, 1/4]';
s.front    = [1/2,3/8, 1/2, 1/2, 3/8, 1/2, 1/4]';

s.volTot   = (s.back + s.front);
s.back     = s.back  ./ s.volTot;
s.front    = s.front ./ s.volTot;
s.upDotN   = s.zx(s.up)  .*s.nx + s.zy(s.up)  .*s.ny  ;
s.downDotN = s.zx(s.down).*s.nx + s.zy(s.down).*s.ny ;
s.nCV      = max([s.from;s.to]);
s.nMC      = length(s.from);
s.nInter   = length(s.nx);
s.nEq      = 2*s.nCV + s.nMC;
s.friction = 0.01;
s.LoD      = 1;

s.theta = atan(s.zy./s.zx);
s.g     = 9.81*cos(s.theta + pi/2);

[s.Ccv,s.Cmc,s.Cinter,s.iInter] = GetSummationMatrices([s.from,s.to],[s.up,s.down],[s.upDotN,s.downDotN]);

% s.Ccv(7,1) = -1;

s.sRho  = 0;
s.sRhov = 0;
s.sRhoe = 10E6*[-1;0;0;1;0;0];

s.dt = 0.20;
rho0 = 9.965569351080000e+02;
e0   = 1.125536123942350e+05;
v0   = 0.01;
s.rho0  = rho0*ones(s.nCV,1);
s.rhoe0 = rho0*e0*ones(s.nCV,1);
s.rhov0 = rho0*v0*ones(s.nMC,1);

s.rhov0(1) = 10*s.rhov0(1);

s.q0 = [s.rho0;s.rhoe0;s.rhov0];

s.iRho  = (1:s.nCV)';
s.iRhoe = s.nCV+ s.iRho ;
s.iRhov = (2*s.nCV + (1:s.nMC))';


% Thermodynamic properites
s.e = s.rhoe0 ./ s.rho0     ;
s.T = Temperature(s.rho0,s.e);
s.P = Pressure(s.rho0,s.T)   ;


s.epsilon = 5E-4;

upwind = IntrepidTwilight.AdamantWave.Semidiscretizations.Quasi2DUpwind(s);
df     = upwind.blockDiagonalJacobian(s.q0);

%{

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
