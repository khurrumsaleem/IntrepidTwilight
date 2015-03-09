function ie = implicitEuler(problem)
    
    qStore = 0                                      ;
    qLast  = problem.initialState.q0                ;
    f      = problem.semidiscretization.closure.rhs ;
    
    ie.deltaq  = @(q,dt) deltaq(q,dt);
    ie.qUpdate = @(qUpdate) updateQ(qUpdate) ;
    ie.qLast   = @() getQLast()  ;
    ie.qStore  = @() getQStore() ;
    
    function dq = deltaq(q,dt)
        dq = dt*f(q);
    end
    function [] = updateQ(q)
        qLast = q;
    end


    function value = getQStore()
        value = qStore;
    end
    function value = getQLast()
        value = qLast;
    end
    
end