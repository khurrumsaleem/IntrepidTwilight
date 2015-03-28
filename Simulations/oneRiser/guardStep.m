function dq = guardStep(q,dq,qLo,qHi)
    
    qGuarded = guardValue(q - dq,qLo,qHi);
    dq       = q - qGuarded;
    
end