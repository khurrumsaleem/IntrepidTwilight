function ie = ImplicitEuler(parameters)
    
%     qStore = 0;
    f      = parameters.f   ;
    dt     = parameters.dt  ;
    qStore = 0              ;
    qLast  = 0              ;
    
    ie.deltaq      = @(q) deltaq(q);
    ie.qLastUpdate = @(qUpdate) updateQLast(qUpdate) ;
    ie.qLast       = @() getQLast()  ;
    ie.qStore      = @() getQStore() ;
    
    function dq = deltaq(q)
        dq = dt*f(q);
    end
    function [] = updateQLast(qNow)
        qLast = qNow;
    end


    function value = getQStore()
        value = qStore;
    end
    function value = getQLast()
        value = qLast;
    end
    
end