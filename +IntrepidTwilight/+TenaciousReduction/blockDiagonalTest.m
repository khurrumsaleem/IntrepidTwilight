clear();
clc();

n = 100;
A = rand(n)         ;
B = rand(n)         ;
C = rand(n+5)       ;

Deco = [[A;B],zeros(2*n,5);C];
D = blkdiag(A,B,C)  ;

xTrue = rand(3*n+5,1)   ;
b     = D*xTrue     ;

tic;
for k = 1:100
    xDirect = D\b;
end
toc;
tic;
for k = 1:100
    xEco = IntrepidTwilight.solvers.blockDiagonalEconomy(Deco,b,[n,n,n+5]);
end
toc;

Show([norm(xTrue),norm(xDirect),norm(xEco)]);