function J = CalculateJacobian(F,x0,epsilons)
    
    Nvars = length(x0);
    dx    = zeros(Nvars,1);

%     dx(1) = epsilons(1);
%     J     = (F(x0 + dx) - F(x0 - dx))/(2*epsilons(1));
%     J     = J(:,ones(1,Nvars));
%     
%     for k = 2:Nvars;
%         dx(k-1) = 0;
%         dx(k)   = epsilons(k);
%         J(:,k)  = (F(x0 + dx) - F(x0 - dx))/(2*epsilons(k));
%     end

    F0    = F(x0);
    dx(1) = epsilons(1);
    J     = (F(x0 + dx) - F0)/(epsilons(1));
    J     = J(:,ones(1,Nvars));
    for k = 2:Nvars;
        dx(k-1) = 0;
        dx(k)   = epsilons(k);
        J(:,k)  = (F(x0 + dx) - F0)/(epsilons(k));
    end
    
end