function r = Residual(timeDiscretization)

    ts                      = timeDiscretization            ;
    r                       = @(q) value(q)                 ;
    r.update                = @(q,t,dt) update(q,t,dt)      ;
    r.blockDiagonalJacobian = @(q) blockDiagonalJacobian(q) ;
    r.jacobian              = @(q) jacobian(q)              ;
    

    %   Residual value
    function r = value(q)
        r = q - ts.qStar();
    end



    %   Update functions
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
    
end