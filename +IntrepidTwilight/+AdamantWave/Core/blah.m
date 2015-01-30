clc();
clear();

B = rand(4);
A = blkdiag(B,B,B,B,B,B,B,B,B,B,B,B);
x = rand(length(A),1);
b = A*x;

tic;
xMat = A\b;
toc;

tic;
blocks    = 12;
blockSize = 4;
xMine     = b;
for k = 1:blocks
    I = ((k-1)*blockSize+1):(k*blockSize);
    xMine(I) = A(I,I)\b(I);
end
toc;

disp(norm(xMat  - x));
disp(norm(xMine - x));

