function x = blockDiagonalEconomy(D,b,blockSize)
    
    x(size(D,1),1) = 0;
    
    %   First block
    s = blockSize(1)    ;
    J = 1:s             ;
    I = J               ;
    x(I) = D(I,J)\b(I)  ;
    
    %   k-th block
    for k = 2:numel(blockSize)
        J    = 1:blockSize(k)   ;
        I    = s + J            ;
        x(I) = D(I,J)\b(I)      ;
        s    = I(end)           ;
    end
    
end