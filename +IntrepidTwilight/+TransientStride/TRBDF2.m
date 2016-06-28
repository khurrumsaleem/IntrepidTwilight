function trbdf2 = TRBDF2(config)

    %   Local properties
    qLast  = 0  ;
    fLast  = 0  ;
    t      = 0  ;
    dt     = 0  ;
    sd     = 0  ;

    c(5,1) = 0  ;
    gamma  = 0  ;
    alpha  = 0  ;


    %   Inherit
    trbdf2 = IntrepidTwilight.executive.Component();
    trbdf2 = trbdf2.changeID(trbdf2,'trbdf2','timediscretization');
    
    %   Define default parameter values
    trbdf2.set('gamma', 2 - sqrt(2));
    trbdf2.set('alpha',      0     );


    %   Imbalance value
    trbdf2.qStar = @(q) qStar(q);
    function value = qStar(q)
        qn1   = q(1:end/2)                ;
        qng   = q(end/2+1:end)            ;
        fng   = sd.rhs(qng , t + gamma*dt);
        fn1   = sd.rhs(qn1 , t +       dt);
        value = [   c(1) * qLast(1) +      c(2) * qng   + dt * c(3) * fn1  ;
                           qLast(1) + dt * c(4) * fLast + dt * c(5) * fng  ];
    end


    %   Late bind
    trbdf2.set('dependencies',{'spacediscretization'});
    trbdf2.bind = @(sd) bind(sd);
    function [] = bind(object)
        if isstruct(object) && object.is('spacediscretization')
            sd = object;
        end
    end
    
    
    %   Prepare
    trbdf2.prepare = @(varargin) prepare(varargin{:});
    function [] = prepare(q0,time,varargin)
        
        %   Cascade preparation
        sd.prepare(q0,time,varargin{:});
        
        %   Initialize
        qLast = q0              ;
        t     = time            ;
        fLast = sd.rhs(qLast,t) ;
        
        %   Set integration weights
        gamma = trbdf2.get('gamma')                                 ;
        alpha = trbdf2.get('alpha')                                 ;
        c(1)  = -alpha*(1-gamma)^2/(gamma* (alpha*(1-gamma) + 1))   ;
        c(2)  = (alpha*(1/gamma-1)+1)/(alpha*(1-gamma)+1)           ;
        c(3)  = (1-gamma)/(alpha*(1-gamma)+1)                       ;
        c(4)  = gamma*alpha/2                                       ;
        c(5)  = gamma*(1-alpha/2)                                   ;

    end

    %   Update function
    trbdf2.update = @(q,t,dt) update(q,t,dt);
    function [] = update(q,time,step)
        %   Update stored states
        qLast  = q(1:end/2)     ;
        fLast  = sd.rhs(qLast,t);
        t      = time           ;
        dt     = step           ;
    end



    %   Accessors
    trbdf2.qLast  = @() getQLast()  ;
    trbdf2.qStore = @() getQStore() ;
    function value = getQStore()
        value = qLast;
    end
    function value = getQLast()
        value = qLast;
    end



    %   Jacobians
    trbdf2.jacobian              = @(q) jacobian(q)             ;
    trbdf2.blockDiagonalJacobian = @(q) blockDiagonalJacobian(q);
    function dq = blockDiagonalJacobian(q)
        dq = sd.blockDiagonalJacobian(q);
    end
    function dq = jacobian(q)
        dq = sd.jacobian(q);
    end
    
    
    if (nargin > 0) && isstruct(config)
        trbdf2.set(config);
    end

end