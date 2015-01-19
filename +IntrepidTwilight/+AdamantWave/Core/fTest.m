clc;
clear();

r2 = cos(pi/4);
s2 = sqrt(2);

s.from     = [1:8, 7,10:11,12,9 ,13:18,19]';
s.to       = [2:9,10,11:12,18,13,14:19, 1]';
s.zx       = [ 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, -1, -1 ]';
s.zy       = [-1,-1,-1,-1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1,  0,  0,  0,  0 ]';
%
s.up       = [1:7               , 6, 7,  8, 9:11,12,12,13:19,20]';
s.down     = [2:8               , 9, 9, 13,10:12,18,19,14:20, 1]';
s.nx       = [ 0, 0, 0, r2,1,1,1,r2,-r2,r2,0,0,0,r2,-r2,0,0,0,-r2,-1,-1,-1,-r2]';
s.ny       = [-1,-1,-1,-r2,0,0,0,r2, r2,r2,1,1,1,r2, r2,1,1,1, r2, 0, 0, 0,-r2]';
%
s.volFrom  = [0.5,0.5,0.5,0.5,0.5,0.5,3/7,0.5,1/3,0.5,0.5,2/3,0.5,0.5,0.5,0.5,0.5,0.5,3/7,0.5]';
s.Ainter   = [1,1,1,s2,1,0.5,1,r2,r2,s2,1,1,1,r2,r2,1,1,1,r2,1,0.5,1,s2]';
s.upDotN   = s.zx(s.up)  .*s.nx + s.zy(s.up)  .*s.ny  ;
s.downDotN = s.zx(s.down).*s.nx + s.zy(s.down).*s.ny ;
s.nCV      = max([s.from;s.to]);
s.nMC      = length(s.from);
s.nInter   = length(s.nx);
s.fFric    = 0.01;
s.LoD      = 1;
s.volTo    = 1 - s.volFrom;

s.theta = atan(s.zy./s.zx);
s.g     = 9.81*cos(s.theta + pi/2);

[s.Ccv,s.Cmc,s.Cp,s.inter] = GetSummationMatrices([s.from,s.to],[s.up,s.down],[s.upDotN,s.downDotN]);

s.Srho  = 0;
s.Srhoe = [0,zeros(1,8),0,0,0,0,zeros(1,6)]';

s.dt   = 0.1;
rho0 = 9.965569351080000e+02;
e0   = 1.125536123942350e+05;
v0   = 0.01;
rho  = rho0*ones(s.nCV,1);
rhoe = rho0*e0*ones(s.nCV,1);
rhov = rho0*v0*ones(s.nMC,1);


s.qOld = [rho;rhoe;rhov];
% f    = @(q) SemidiscreteUpwind(q,s);

% Thermodynamic properites
s.rho = rho;
s.e = rhoe ./ s.rho       ;
s.T = Temperature(s.rho,s.e);
s.P = Pressure(s.rho,s.T)   ;


q = SolveClosedLoopSystem(s.qOld,s);

% while true
%     
%     r    = @(q) (q-qOld)-dt*f(q);
%     qNew = JFNKHouseholder(1.0008*qOld,r,1E-8);
%     
% end





