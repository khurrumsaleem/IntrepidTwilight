function xMin = goldenSectionSearch(f,a,b,tolerance)
    
    narginchk(3,4);
    
    if (nargin < 3) || isempty(tolerance)
        tolerance = 1E-8*max([a,b]);
    end
    
    %   Interval set-up
    phi = (sqrt(5) - 1)/2   ;
    c   = b + phi*(a-b)     ;
    d   = a - phi*(a-b)     ;
    fc  = f(c)              ;
    fd  = f(d)              ;
    
    while (abs(c - d) > tolerance)
        
        if fc < fd
            
            %   Shift interval bounds
            b  = d;     %fb = fd;
            d  = c;     fd = fc;
            
            %   Sample new point
            c  = b + phi*(a-b)  ;
            fc = f(c)           ;

        else
            
            %   Shift interval bounds
            a  = c;     %fa = fc;
            c  = d;     fc = fd;
            
            %   Sample new point
            d  = a - phi*(a-b)  ;
            fd = f(d)           ;

        end

    end
    
    xMin = c;
    
end