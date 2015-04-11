function dq = guardStep(q,dq,qLo,qHi)

    if any(abs(dq) > 1E6)
        g = [];
    end
    qGuarded = guardValue(q - dq,qLo,qHi);
    dq       = q - qGuarded;
    
end