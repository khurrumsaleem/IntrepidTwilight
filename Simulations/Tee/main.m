clc();
clear();



T1  = 300                       ;   %   [K]
P1  = 120E3                     ;   %   [Pa]
P7  = 119E3                     ;   %   [Pa]
P10 = 120E3 - 9.81*996.56*0.3   ;   %   [Pa]
%
%   Setup
hem = testHEM(T1,P1,P7,P10);
hem.set('evolver','time.span'         , [0,5000] );
hem.set('evolver','time.step.maximum' ,    1    );
hem.set('evolver','time.step.minimum' ,   1E-12  );
hem.set('evolver','time.step.goal'    ,   1     );
hem.set('evolver','saveRate'          ,   0.5    );
%
%   Run
trun = tic;
hem.run();
tElapsed(k) = toc(trun);
[q{k},Dq{k},t{k}] = hem.results()       ;