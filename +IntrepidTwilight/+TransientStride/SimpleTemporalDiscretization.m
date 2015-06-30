function ts = SimpleTemporalDiscretization(spatialDiscretization)

    %   Local properties
    qLast  = 0                      ;
    dt     = 0                      ;
    sd     = spatialDiscretization  ;

    %   Methods
    ts.qStar                 = @(q) qStar(q)                ;
    ts.qUpdate               = @(q,t,dt) update(q,t,dt)     ;
    ts.qLast                 = @() getQLast()               ;
    ts.qStore                = @() getQStore()              ;
    ts.jacobian              = @(q) jacobian(q)             ;
    ts.blockDiagonalJacobian = @(q) blockDiagonalJacobian(q);


    %   Imbalance function
    function value = qStar(q)
        value = qLast + dt * sd.rhs(q);
    end



    %   Update function
    function [] = update(q,t,step)
        %   Update stored states
        qLast  = q      ;
        dt     = step   ;
        
        %   Update spatial discretization's time
        sd.update(t,step);
    end



    %   Accessors
    function value = getQStore()
        value = qLast;
    end
    function value = getQLast()
        value = qLast;
    end



    %   Jacobians
    function dq = blockDiagonalJacobian(q)
        dq = sd.blockDiagonalJacobian(q);
    end
    function dq = jacobian(q)
        dq = sd.jacobian(q);
    end

end