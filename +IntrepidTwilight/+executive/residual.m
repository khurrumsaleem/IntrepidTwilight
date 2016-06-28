function r = Residual(config)

    %   Inherit and setup
    r = IntrepidTwilight.executive.Component();
    r = r.changeID(r,'residual','residual');


    %   Binder
    ts     = [];
    r.set('dependencies',{'timediscretization'});
    r.bind = @(object) bind(object);
    function [] = bind(object)
        if isstruct(object) && object.is('timediscretization')
            ts = object;
        end
    end


    %   Expansion point
    expansionPoint   = [];
    r.expansionPoint = @() getExpansionPoint();
    function q = getExpansionPoint()
        q = expansionPoint;
    end


    %   Residual value
    r.value = @(q) value(q);
    function r = value(q)
        expansionPoint = q;
        r = q - ts.qStar(q);
    end


    %   Updater
    r.update = @(q,t,dt) update(q,t,dt);
    function [] = update(q,t,dt)
        ts.update(q,t,dt);
    end


    %   Jacobians
    r.blockDiagonalJacobian = @(q) blockDiagonalJacobian(q) ;
    r.jacobian              = @(q) jacobian(q)              ;
    function drdq = blockDiagonalJacobian(q)
        dfdq = ts.blockDiagonalJacobian(q)                              ;
        drdq = cellfun(@(c) eye(size(c)) - c,dfdq,'UniformOutput',false);
    end
    function drdq = jacobian(q)
        dfdq = ts.jacobian(q)           ;
        drdq = eye(size(dfdq)) - dfdq   ;
    end
    
    
    
    
    
    %
    %   Prepare:
    %       Bind step and value functions
    r.prepare      = @(varargin) prepare(varargin{:})   ;
    guardStep_     = []                                 ;
    guardValue_    = []                                 ;
    isNotPreapared = true                               ;
    function [] = prepare(q,t,varargin)
        if isNotPreapared
        %   Cascade preparation
        ts.prepare(q,t,varargin{:});
        %
        expansionPoint = ts.qLast();
        %
        guardStep_  = r.get('guard.step');
        guardValue_ = r.get('guard.value');
        %
        isNotPreapared = false;
        end
    end




    %   Default guards do nothing
    r.guard.step  = @(q,dq) guardStep(q,dq) ;
    r.guard.value = @(q)    guardValue(q)   ;
    r.set('guard.value',@(q)     q)  ;
    r.set('guard.step' ,@(q,dq) dq)  ;
    %
    %   Wraps r in the closure which is passed to the, possibly user-defined, guards
    function dq = guardStep(q,dq)
        dq = guardStep_(q,dq);
    end
    function q = guardValue(q)
        q     = guardValue_(q);
    end




    %   Set key-store 
    if (nargin >= 1)
        r.set(config);
    end

end