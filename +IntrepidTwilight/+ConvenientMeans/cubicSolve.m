function r = cubicSolve(p)
    
    c2 = p(:,2)./p(:,1);
    c1 = p(:,3)./p(:,1);
    c0 = p(:,4)./p(:,1);
    
    Q = (3*c1 - c2.^2)/9                ;
    R = (9*c2.*c1 - 27*c0 - 2*c2.^3)/54 ;
    D = sqrt(Q.^3 + R.^2)               ;
    S = sign(R+D).*abs(R + D).^(1/3)    ;
    T = sign(R-D).*abs(R - D).^(1/3)    ;

    SpT = S+T               ;
    SmT = 1i*sqrt(3)*(S-T)/2;
    r   = bsxfun(@plus,-c2/3,[SpT,-SpT/2+SmT,-SpT/2-SmT]);
    
end