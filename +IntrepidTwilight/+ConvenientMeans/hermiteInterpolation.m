function y = hermiteInterpolation(xs,ys,Dys,x)
    
    %   Filter out NaNs which indicates data desired not to be used
    iC0   = not(isnan(ys(1,:)))     ;
    iC1   = not(isnan(Dys(1,:)))    ;
    xC0   = xs(iC0)                 ;
    xC1   = xs(iC1)                 ;
    ys    = ys(:,iC0)               ;
    Dys   = Dys(:,iC1)              ;
    [m,~] = size(ys)                ;
    n     = nnz(iC0) + nnz(iC1) - 1 ;

    %   Scale
    mu    = mean(xs);
    sigma = std(xs);
    xC0   = (xC0 - mu)/sigma;
    xC1   = (xC1 - mu)/sigma;
    
    %   Build the block Vandermonde matrix
    V  =[   bsxfun(@power,xC0,n:-1:0);
            [bsxfun(@times,bsxfun(@power,xC1,n-1:-1:1),n:-1:2),xC1*0+1,xC1*0]/sigma...
        ].';
    % alpha = [ones(1,npoly+1);zeros(nherm,npoly+1)];
    % if (nherm>=1)
    %     alpha(2,:) = 0:npoly;
    %     for k = 3:nherm+1
    %         alpha(k,:) = alpha(k-1,:) .* (alpha(2,:) - k+2);
    %     end
    % end


    %   Solve for coefficients
    c(m,n+1) = 0;
    for k = 1:m
        c(k,:) = [ys(k,:),Dys(k,:)] / V;
    end


    %   Evaluate Polynomial at requested point using Horner's Method
    y = c(:,1)              ;
    x = (x(:).' - mu)/sigma ;
    for k = 2:n+1
        y = bsxfun(@plus,c(:,k),bsxfun(@times,y,x));
    end

end


