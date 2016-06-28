function r = cubicSolve2(p)
    
    %   Reduce cubic
    p   = bsxfun(@rdivide,p(:,2:4),p(:,1))          ;
    tau = p(:,1).^2/3 - p(:,2)                      ;
    ksi = 2/27*p(:,1).^3 - p(:,1).*p(:,2)/3 + p(:,3);
    
    %   Discriminant tests
    delta        = 4*tau.^3 - 27*ksi.^2 ;
    threeAreReal = delta >= 0           ;
    oneIsReal    = not(threeAreReal)    ;
    
    %   Allocate
    r = p;
    
    %   Calculate roots
    %
    %   Case 1
    mask      = threeAreReal;
    if any(mask)
        r(mask,:) = calculateThreeReal(tau(mask),ksi(mask),delta(mask));
        r(mask,:) = bsxfun(@plus,-p(mask,1)/3,r(mask,:));
    end
    %
    %   Case 2
    mask      = oneIsReal & tau < 0;
    if any(mask)
        r(mask,:) = calculateOneRealNegative(tau(mask),ksi(mask));
        r(mask,:) = bsxfun(@plus,-p(mask,1)/3,r(mask,:));
    end
    %
    %   Case 3
    mask      = oneIsReal & tau > 0;
    if any(mask)
        r(mask,:) = calculateOneRealPositive(tau(mask),ksi(mask));
        r(mask,:) = bsxfun(@plus,-p(mask,1)/3,r(mask,:));
    end
    
    
end

function r = calculateThreeReal(tau,ksi,delta)
    
    rho   = sqrt(4/3*tau);
    omega = atan(9*ksi./sqrt(3*delta))/3;
    r     = [rho.*cos(omega+pi/6),-rho.*cos(omega-pi/6),rho.*sin(omega)];
    
end

function r = calculateOneRealNegative(tau,ksi)
    
    rho   = 1i*sqrt(4/3*tau);
    phi   = asinh(-4*ksi./rho.^3)/3;
    omega = 1i*sqrt(3)*cosh(phi);
    r     = [rho.*sinh(phi),-rho.*(omega+sinh(phi))/2,rho.*(omega-sinh(phi))/2];
    
end

function r = calculateOneRealPositive(tau,ksi)
    
    rho   = sqrt(4/3*tau);
    phi   = acosh(-4*ksi./rho.^3)/3;
    omega = 1i*sqrt(3)*sinh(phi);
    r     = [rho.*cosh(phi),-rho.*(omega+cosh(phi))/2,rho.*(omega-cosh(phi))/2];
    
end


