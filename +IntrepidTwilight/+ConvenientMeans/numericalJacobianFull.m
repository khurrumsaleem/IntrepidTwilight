function J = numericalJacobianFull(f,x0,epsilons)
    
    
    if (nargin < 3) || isempty(epsilons)
        aboveOne = abs(x0) >= 1;
        epsilons = 1E-7 * (x0.*aboveOne + not(aboveOne));
    end
    
    Nvars       = length(x0);
    dx(Nvars,1) = 0         ;

    f0    = f(x0);
    dx(1) = epsilons(1);
    J     = (f(x0 + dx) - f0)/(epsilons(1));
    J     = J(:,ones(1,Nvars));
    for k = 2:Nvars;
        dx(k-1) = 0;
        dx(k)   = epsilons(k);
        J(:,k)  = (f(x0 + dx) - f0)/(epsilons(k));
    end
    
end