clc;
clear();
import IntrepidTwilight.AdamantWave.Semidiscretizations.Quasi2DUpwind.*;

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
s.volFrom  = [1/2,1/2, 3/8, 1/2, 1/2, 3/8, 1/4]';
s.volTo    = [1/2,3/8, 1/2, 1/2, 3/8, 1/2, 1/4]';

s.volTot   = (s.volFrom + s.volTo);
s.volFrom  = s.volFrom./s.volTot;
s.volTo    = s.volTo  ./s.volTot;
s.upDotN   = s.zx(s.up)  .*s.nx + s.zy(s.up)  .*s.ny  ;
s.downDotN = s.zx(s.down).*s.nx + s.zy(s.down).*s.ny ;
s.nCV      = max([s.from;s.to]);
s.nMC      = length(s.from);
s.nInter   = length(s.nx);
s.nEq      = 2*s.nCV + s.nMC;
s.fFric    = 0.01;
s.LoD      = 1;

s.theta = atan(s.zy./s.zx);
s.g     = 9.81*cos(s.theta + pi/2);

[s.Ccv,s.Cmc,s.Cp,s.inter] = GetSummationMatrices([s.from,s.to],[s.up,s.down],[s.upDotN,s.downDotN]);

% s.Ccv(7,1) = -1;

s.Srho  = 0;
s.Srhoe = 10E6*[-1;0;0;1;0;0];

s.dt = 0.20;
rho0 = 9.965569351080000e+02;
e0   = 1.125536123942350e+05;
v0   = 0.01;
rho  = rho0*ones(s.nCV,1);
rhoe = rho0*e0*ones(s.nCV,1);
rhov = rho0*v0*ones(s.nMC,1);

rhov(1) = 10*rhov(1);

s.rhoMask  = (1:s.nCV)';
s.rhoeMask = s.nCV+ s.rhoMask ;
s.rhovMask = (2*s.nCV + (1:s.nMC))';

s.internalMask = [s.rhoMask(1:end);s.rhoeMask(1:end);s.rhovMask(1:end)];

q0      = [rho;rhoe;rhov]   ;
qstar   = q0                ;
s.qstar = qstar             ;
% f    = @(q) SemidiscreteUpwind(q,s);

% Thermodynamic properites
s.rho = rho;
s.e = rhoe ./ s.rho         ;
s.T = Temperature(s.rho,s.e);
s.P = Pressure(s.rho,s.T)   ;

% r    = @(qND) (qND - q0./qstar) - s.dt*SemidiscreteUpwind(qND.*qstar,s)./qstar;
r    = @(q) (q - q0) - s.dt*SemidiscreteUpwind(q,s);
drdq = CalculateJacobian(r,q0,1E-8*q0);


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
