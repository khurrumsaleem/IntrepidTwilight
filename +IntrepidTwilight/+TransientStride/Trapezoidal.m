function ts = Trapezoidal(config)

    %   Local properties
    qLast = 0   ;
    fLast = 0   ;
    dt    = 0   ;
    sd    = 0   ;

    %   Inherit
    ts = IntrepidTwilight.executive.Component();
    ts = ts.changeID(ts,'trapezoidal','timediscretization');


    %   Imbalance value
    ts.qStar = @(q) qStar(q);
    function value = qStar(q)
        value = qLast + dt/2 * (fLast + sd.rhs(q));
    end


    %   Late bind
    ts.set('dependencies',{'spacediscretization'});
    ts.bind = @(sd) bind(sd);
    function [] = bind(object)
        if isstruct(object) && object.is('spacediscretization')
            sd = object;
        end
    end
    
    
    %   Prepare
    ts.prepare = @() prepare();
    function [] = prepare()
        %   Cascade preparation
        sd.prepare();
        %
        qLast = sd.getAll();
        %
    end



    %   Update function
    ts.update = @(q,t,dt) update(q,t,dt);
    function [] = update(q,t,step)
        %   Update stored states
        qLast  = q              ;
        fLast  = sd.rhs(qLast)  ;
        dt     = step           ;
        
        %   Update spatial discretization's time
        sd.update(t);
    end



    %   Accessors
    ts.qLast  = @() getQLast()  ;
    ts.qStore = @() getQStore() ;
    function value = getQStore()
        value = qLast;
    end
    function value = getQLast()
        value = qLast;
    end



    %   Jacobians
    ts.jacobian              = @(q) jacobian(q)             ;
    ts.blockDiagonalJacobian = @(q) blockDiagonalJacobian(q);
    function dq = blockDiagonalJacobian(q)
        dq = sd.blockDiagonalJacobian(q);
    end
    function dq = jacobian(q)
        dq = sd.jacobian(q);
    end
    
    
    if (nargin > 0) && isstruct(config)
        ts.set(config);
    end

end