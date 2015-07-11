function ts = SimpleTemporalDiscretization(spatialDiscretization)

    %   Local properties
    qLast = 0   ;
    dt    = 0   ;
    sd    = 0   ;
    

    %   Bind at construction if passed
    if (nargin >= 1)
        sd = spatialDiscretization ;
    end

    %   Methods
    ts.is                    = @(s) strcmpi(s,'timediscretization') ;
    ts.set                   = @(type,object) set(type,object)      ;
    ts.qStar                 = @(q) qStar(q)                        ;
    ts.update                = @(q,t,dt) update(q,t,dt)             ;
    ts.qLast                 = @() getQLast()                       ;
    ts.qStore                = @() getQStore()                      ;
    ts.jacobian              = @(q) jacobian(q)                     ;
    ts.blockDiagonalJacobian = @(q) blockDiagonalJacobian(q)        ;


    %   Imbalance value
    function value = qStar(q)
        value = qLast + dt * sd.rhs(q);
    end



    %   Late bind
    function [] = set(type,object)
        switch(lower(type))
            case('spacediscretization')
                if isstruct(object) && object.is('spacediscretization')
                    sd = object;
                end
        end
    end



    %   Update function
    function [] = update(q,t,step)
        %   Update stored states
        qLast  = q      ;
        dt     = step   ;
        
        %   Update spatial discretization's time
        sd.update(t);
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