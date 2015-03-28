clc();
clear();

n = 10;
x = rand(n,1);
A = rand(n,n);
b = A*x;
r = @(x) exp(x) - (1:n)';
% r = @(x) x.^10;
drdq = @(x) diag(exp(x));
x0 = ones(n,1);

% x   = x0;
% err = r(x);
% while norm(err,2) > 1E-12
%     dx  = drdq(x)\err   ;
%     x   = x - dx        ;
%     err = r(x)          ;
% end


problem.solver.name = 'JFNK';
solver = IntrepidTwilight.executive.Solver(problem);
solver.gmres.nu = 0.85;

N = 1;
tic;
for k = 1:N
    [xsol1,stats1] = IntrepidTwilight.TenaciousReduction.JFNKnew1(x0,r,solver);
end
toc;
tic;
for k = 1:N
[xsol15,stats15] = IntrepidTwilight.TenaciousReduction.JFNKnew15(x0,r,solver);
end
toc;
tic;
for k = 1:N
[xsol3,stats3] = IntrepidTwilight.TenaciousReduction.JFNKnew3(x0,r,solver);
end
toc;
tic;
for k = 1:N
[xsol4,stats4] = IntrepidTwilight.TenaciousReduction.JFNKnew4(x0,r,solver);
end
toc;

Show([xsol1(1);xsol15(1);xsol3(1);xsol4(1)]);

semilogy(...
    (1:stats1.iterations)',stats1.norm,'-*',...
    (1:stats15.iterations)',stats15.norm,'-o',...
    (1:stats3.iterations)',stats3.norm,'--o',...
    (1:stats4.iterations)',stats4.norm,'--o');
legend('1','1.5','3','4');
