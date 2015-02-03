clc();
clear();

B = rand(4);
A = blkdiag(B,B,B,B,B,B,B,B,B,B,B,B);
D = [B;B;B;B;B;B;B;B;B;B;B;B];
x = rand(length(A),1);
b = A*x;

for k = 1:1000
    xMat  = mldivide(A,b);
    xMine = solveSquareBlockSystem(D,b);
end

Show([norm(xMat  - x),norm(xMine - x)]');

