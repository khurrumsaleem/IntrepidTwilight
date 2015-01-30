clc();
clear();

N = 100;
r = @(x) exp(x) - (1:N)';

x = JFNKHouseholder(zeros(N,1),r,1E-8);

AbsError = abs(log(1:N)'-x) + eps();
semilogy(AbsError);