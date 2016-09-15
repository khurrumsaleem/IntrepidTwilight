function y = linearlyImplicitEuler(rhs,y0,t)
    
    nt      = numel(t)          ;
    y       = repmat(y0,1,nt)   ;
    epsilon = 1E-7              ;
    
    for k = 2:nt
        %
        yk   = y(:,k-1)         ;
        dt   = t(k) - t(k-1)    ;
        rhsk = rhs(yk,t(k-1))   ;
        %
        b      = dt*rhsk                                                    ;
        Afun   = @(dy) dy - dt * (rhs(yk+epsilon*dy,t(k-1)) - rhsk)/epsilon ;
        dy     = IntrepidTwilight.ConvenientMeans.GMRES(Afun,b,y0*0,@(x)x)  ;
        y(:,k) = y(:,k-1) + dy                                              ;
    end

end