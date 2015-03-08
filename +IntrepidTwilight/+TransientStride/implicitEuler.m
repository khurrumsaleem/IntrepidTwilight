function ie = implicitEuler(f,dt)
    
    qStore = 0          ;
    qLast  = 0          ;
    
    ie.deltaq  = @(q) deltaq(q);
    ie.qUpdate = @(qUpdate) updateQ(qUpdate) ;
    ie.qLast   = @() getQLast()  ;
    ie.qStore  = @() getQStore() ;
    
    function dq = deltaq(q)
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