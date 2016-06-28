function rk4 = RungeKutta4(config)

    %   Local properties
    qLast = 0   ;
    time  = 0   ;
    dt    = 0   ;
    sd    = 0   ;

    %   Inherit
    rk4 = IntrepidTwilight.executive.Component();
    rk4 = rk4.changeID(rk4,'RungeKutta4','timediscretization');


    %   Imbalance value
    rk4.qStar = @(q) qStar(q);
    function value = qStar(q)
        
        dt = 0.5 * dt;
        
        k1 = sd.rhs(q);
        sd.update(time + dt/2);
        k2 = sd.rhs(q + dt/2*k1);
        k3 = sd.rhs(q + dt/2*k2);
        sd.update(time + dt);
        k4 = sd.rhs(q + dt*k3);
        
        value = q + dt/6*(k1 + 2*k2 + 2*k3 + k4);
    end


    %   Late bind
    rk4.set('dependencies',{'spacediscretization'});
    rk4.bind = @(sd) bind(sd);
    function [] = bind(object)
        if isstruct(object) && object.is('spacediscretization')
            sd = object;
        end
    end
    
    
    %   Prepare
    rk4.prepare = @() prepare();
    function [] = prepare()
        %   Cascade preparation
        sd.prepare();
        %
        qLast = sd.getAll();
        %
    end



    %   Update function
    rk4.update = @(q,t,dt) update(q,t,dt);
    function [] = update(q,t,step)
        %   Update stored states
        qLast = q       ;
        time  = t       ;
        dt    = step    ;
        
        %   Update spatial discretization's time
        sd.update(t);
    end



    %   Accessors
    rk4.qLast  = @() getQLast()  ;
    rk4.qStore = @() getQStore() ;
    function value = getQStore()
        value = qLast;
    end
    function value = getQLast()
        value = qLast;
    end



    %   Jacobians
    rk4.jacobian              = @(q) jacobian(q)             ;
    rk4.blockDiagonalJacobian = @(q) blockDiagonalJacobian(q);
    function dq = blockDiagonalJacobian(q)
        dq = sd.blockDiagonalJacobian(q);
    end
    function dq = jacobian(q)
        dq = sd.jacobian(q);
    end
    
    
    if (nargin > 0) && isstruct(config)
        rk4.set(config);
    end

end