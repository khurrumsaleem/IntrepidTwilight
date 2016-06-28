function y = hermiteInterpolation2(xs,ys,Dys,x)
    
    [m,n] = size(ys);
    n     = 2*n - 1;

    %   Build the block Vandermonde matrix
    xs = xs(:)  ;
    V  =[   bsxfun(@power,xs,n:-1:0);
            [bsxfun(@mtimes,bsxfun(@power,xs,n-1:-1:1),n:-1:2),xs*0+1,xs*0]...
        ].';


    %   Solve for coefficients
    c(m,n+1) = 0;
    for k = 1:m
        c(k,:) = [ys(k,:),Dys(k,:)] / V;
    end


    %   Evaluate Polynomial at requested point using Horner's Method
    y = c(:,1)              ;
    x = x(:).'              ;
    for k = 2:n+1
        y = bsxfun(@plus,c(:,k),bsxfun(@times,y,x));
    end

end


