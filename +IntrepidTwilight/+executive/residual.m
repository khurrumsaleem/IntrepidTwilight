function r = Residual(timeDiscretization)

    ts = 0;

    %   Bind at construction if passed
    if (nargin >= 1)
        bind('timediscretization',timeDiscretization);
    end


    %   Default guard (does nothing)
    r.guard.step  = @(q,dq) guardStep(q,dq) ;
    r.guard.state = @(q)    guardState(q)   ;


    %   Initialization for user-defined guard
    r.userDefined.guard.step  = [];
    r.userDefined.guard.state = [];


    %   Public methods
    r.type                  = 'residual'                    ;
    r.is                    = @(s) strcmpi(s,r.type)        ;
    r.bind                  = @(object) bind(object)        ;
    r.value                 = @(q) value(q)                 ;
    r.update                = @(q,t,dt) update(q,t,dt)      ;
    r.blockDiagonalJacobian = @(q) blockDiagonalJacobian(q) ;
    r.jacobian              = @(q) jacobian(q)              ;


    %   Residual value
    function r = value(q)
        r = q - ts.qStar(q);
    end



    %   Late binder
    function [] = bind(object)
        if isstruct(object)
            switch(object(1).type)
                case('timediscretization')
                        ts = object;
            end
        end
    end
    


    %   Update function
    function [] = update(q,t,dt)
        ts.update(q,t,dt);
    end



    %   Jacobians
    function drdq = blockDiagonalJacobian(q)
        dfdq = ts.blockDiagonalJacobian(q)                              ;
        drdq = cellfun(@(c) eye(size(c)) - c,dfdq,'UniformOutput',false);
    end
    function drdq = jacobian(q)
        dfdq = ts.jacobian(q)           ;
        drdq = eye(size(dfdq)) - dfdq   ;
    end



    %   Guard
    function [value,dq] = guardStep(q,dq)
        if isa(r.userDefined.guard.step,'function_handle');
            [value,dq] = r.userDefined.guard.step(r,q,dq);
        else
            value = r.value(q-dq);
        end
    end
    function [value,q] = guardState(q)
        if isa(r.userDefined.guard.state,'function_handle');
            [value,q] = r.userDefined.guard.state(r,q);
        else
            value = r.value(q);
        end
    end

end