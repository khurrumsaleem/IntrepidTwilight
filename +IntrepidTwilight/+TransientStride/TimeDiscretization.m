function ts = TimeDiscretization(config)
    
    %   Local properties
    qLast = 0   ;
    t     = 0   ;
    dt    = 0   ;
    sd    = 0   ;
    
    %   Inherit
    ts = IntrepidTwilight.executive.Component();
    ts = ts.changeID(ts,'timediscretization','timediscretization');
    
    
    %   Imbalance value
    ts.qStar = @(q) qStar(q);
    function value = qStar(q)
        value = qLast + dt * sd.rhs(q,t);
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
    ts.prepare    = @(varargin) prepare(varargin{:});
    isNotPrepared = true;
    function [] = prepare(q0,t,varargin)
        if isNotPrepared
            %   Cascade preparation
            sd.prepare(q0,t);
            qLast = sd.getAll();
            
            isNotPrepared = false;
        end
    end
    
    
    
    %   Update function
    ts.update = @(q,t,dt) update(q,t,dt);
    function [] = update(q,time,step)
        %   Update stored states
        qLast  = q      ;
        t      = time   ;
        dt     = step   ;
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
        dq = cellfun(@(c) dt * c,sd.blockDiagonalJacobian(q),'UniformOutput',false);
    end
    function dq = jacobian(q)
        dq = dt * sd.jacobian(q);
    end
    
    
    if (nargin > 0) && isstruct(config)
        ts.set(config);
    end
    
end