function q = guardValue(q,qLo,qHi)
    
    maskLo = q <= qLo;
    maskHi = q >= qHi;
    
    q(maskLo) = 1.0001*qLo(maskLo);
    q(maskHi) = 0.9999*qHi(maskHi);
    
end