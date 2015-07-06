function [pointMin,valueMin] = goldenSectionMinimizeLazily(f,interval,values,tolerance)
    
    narginchk(3,4);
    if (nargin < 3) || isempty(tolerance)
        tolerance = 1E-8*(interval(2)-interval(1));
    end
    
    
    %   Interval set-up
    phi     = (sqrt(5) - 1)/2                                           ;
    a       = interval(1)                                               ;
    b       = interval(2)                                               ;
    points  = [    a     , b + phi*(a-b) , a - phi*(a-b) ,      b   ]   ;
    values  = [values(1) , f(points(2))  , f(points(3))  , values(2)]   ;



    %   If both interior samples are NOT less than the minimum of the end points,
    %   no search is performed.  Either the minimum of the function is the minimum
    %   bounding point since it is Monotonic in some fashion, or the function has
    %   a complex or pathological convexity that will not be explored (lazy).
    if all(values([2,3]) < min(values([1,4])))
        
        %   Iterate
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
        
    end
    
    %   Return minimum value and associated point
    valueMin = min(values)                  ;
    pointMin = points(values == valueMin)   ;
    
end