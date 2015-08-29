function r = Residual(config)

    %   Inherit and setup
    r = IntrepidTwilight.executive.Component();
    r = r.changeID(r,'residual','residual');
    r.dependencies = {'timediscretization'};


    %   Binder
    ts     = [];
    r.bind = @(object) bind(object);
    function [] = bind(object)
        if isstruct(object) && object.is('timediscretization')
            ts = object;
        end
    end


    %   Residual value
    r.value = @(q) value(q);
    function r = value(q)
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


    %   Default guard
    r.guard.step  = @(q,dq) guardStep(q,dq) ;
    r.guard.state = @(q)    guardState(q)   ;
    r.set('guard.state',@(r,q)   simpleBacktracker(r,q));
    r.set('guard.step' ,@(r,q,dq)simpleBacktracker(r,q,dq));
    %
    %   Wraps r in the closure which is passed to the, possibly user-defined, guards
    function [dq,rValue] = guardStep(q,dq)       
        guard = r.get('guard.step');
        [dq,rValue] = guard(r,q,dq);
    end
    function [q,rValue] = guardState(q)
        guard = r.get('guard.state');
        [q,rValue] = guard(r,q);
    end
    %
    %   Only looks for NaNs and relaxs the value
    function [qValue,rValue] = simpleBacktracker(r,q,dq)
        switch(nargin)
            case(2)
                %   State
                rValue = r(q);
                while any(isnan(rValue))
                    q      = 0.9*q  ;
                    rValue = r(q)   ;
                end
                qValue = q;
            case(3)
                    %   Step
                rValue = r(q-dq);
                while any(isnan(rValue))
                    dq     = 0.5*dq     ;
                    rValue = r(q - dq)  ;
                end
                qValue = dq;
        end
    end



    %   Set key-store 
    if (nargin >= 1)
        r.set(config);
    end

end