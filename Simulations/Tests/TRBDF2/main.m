clc();
clear();

syms y(t)

f    = diff(tanh(6*t-10),t);
yana = matlabFunction(dsolve(diff(y,t) == f,y(0) == 0));

a  = 6;
b  = 10;
y0 = 0                  ;
f  = matlabFunction(f)  ;

sd = IntrepidTwilight.AdamantWave.SpaceDiscretization();
sd.rhs    = @(q,t) f(t);

alphas = [0;0.5;1];
n      = numel(alphas);
q      = cell(1,3);
t      = cell(1,3);
tspan  = [0,10];
dt     = diff(tspan)/50;
for k = 1:n
    evolver = makeEvolver(sd,alphas(k));
    evolver.set('time.span'         ,  tspan );
    evolver.set('time.step.maximum' ,    1   );
    evolver.set('time.step.minimum' ,   1E-7 );
    evolver.set('time.step.goal'    ,   dt   );
    evolver.set('initialCondition'  ,   y0   );
    evolver.set('saveRate'          ,   dt   );
    evolver.run();
    [q{k},t{k}] = evolver.getData();
end

tana = linspace(tspan(1),tspan(2),1E3);
qana = yana(tana)       ;

args{1,3*n} = [];
args(1:3:end) = t;
args(2:3:end) = q;
args(3:3:end) = {'o','s','*'};
plot(tana,qana,args{:});
legend('Analytical','\alpha = 0','\alpha = 1/2','\alpha = 1');

