% clc();
% clear();
% 
% N = 1E3;
% A = sprand(N,N,100/N^2);
% T = spdiags([ones(N,1),2*ones(N,1),ones(N,1)],[-1,0,1],N,N);
% A = A + T;
% 
% x = rand(N,1);
% b = A*x;

% [xMat1,~,~,~,rezMat1] = gmres(A,b,50,1E-8,100);
% [xMine1,rezMine1] = GMRESHouseholderCore(A,b,zeros(N,1),50,50*100,8E-5,0.85,@(v)v,@(v)v);

tic;
[xMat2,~,~,~,rezMat2] = gmres(A,b,50,1E-8,100,T);
toc;
tic;
[xMine2,rezMine2] = GMRESHouseholderCore(A,b,zeros(N,1),50,50*100,1.503E-5,0.15,@(v)T\v,@(v)v);
toc;


% figure(1)
%     loglog(1:length(rezMat1),rezMat1,1:length(rezMine1),rezMine1,'o');
figure(2);
    loglog(1:length(rezMat2),rezMat2,1:length(rezMine2),rezMine2,'o');
    
disp([rezMat2(end),rezMine2(end)]);