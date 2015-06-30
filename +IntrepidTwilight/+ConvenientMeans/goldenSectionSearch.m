function [pointMin,valueMin] = goldenSectionSearch(f,interval,values,tolerance)
    
    narginchk(3,4);
    
    if (nargin < 3) || isempty(tolerance)
        tolerance = 1E-8*max(interval);
    end
    
    %   Interval set-up
    phi     = (sqrt(5) - 1)/2                                           ;
    a       = interval(1)                                               ;
    b       = interval(2)                                               ;
    points  = [    a     , b + phi*(a-b) , a - phi*(a-b) ,      b   ]   ;
    values  = [values(1) , f(points(2))  , f(points(3))  , values(2)]   ;
    
    while (abs(points(2)-points(3)) > tolerance)
        
        if values(2) < values(3)
            
            %   Shift interval bounds
            points(4)  = points(3);     values(4) = values(3);
            points(3)  = points(2);     values(3) = values(2);
            
            %   Sample new point
            points(2) = points(4) + phi*(points(1)-points(4))   ;
            values(2) = f(points(2))                            ;

        else
            
            %   Shift interval bounds
            points(1)  = points(2);     values(1) = values(2);
            points(2)  = points(3);     values(2) = values(3);
            
            %   Sample new point
            points(3) = points(1) - phi*(points(1)-points(4))   ;
            values(3) = f(points(3))                            ;

        end

    end
    
    %   Return minimum value and associated point
    valueMin = min(values)                  ;
    pointMin = points(values == valueMin)   ;
    
end