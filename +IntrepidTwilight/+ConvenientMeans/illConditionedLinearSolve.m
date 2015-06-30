function [x,stats] = illConditionedLinearSolve(A,b,x,gamma)

    %   Minor set-up
    At        = A'      ;
    AAt       = A*At    ;
    gamma     = 1-gamma ;
    iter      = 0       ;
    iterMax   = 1000    ;
    tolerance = 1E-13   ;
    r         = A*x - b ;
    
    %   Iterate
    while (norm(r,2) > tolerance) && (iter < iterMax)
        
        %   Update x
        v1    = AAt*r;
        v2    = A*r;
        alpha = (((v1'*r)*v2 - (v2'*r)*v1)'*v1)/(((v2'*r)*v1 - (v1'*r)*v2)'*v2);
        u     = alpha*r + At*r;
        v     = v1 + alpha*v2;
        dx    = gamma * (r'*v)/(v'*v) * u;
        x     = x - dx;
        
        %   Update residual
        r     = A*x - b;

        %   Increment
        iter = iter + 1;
    end
    
    %   Run data struct
    stats.iterations = iter;
    stats.converged  = (norm(r,2) <= tolerance) && (iter < iterMax);
end