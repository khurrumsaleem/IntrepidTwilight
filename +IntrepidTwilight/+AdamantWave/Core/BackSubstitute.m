function x = BackSubstitute(R,b)
    
    N    = length(b)    ;
    x    = b            ;
    x(N) = b(N)/R(N,N)  ;
    
    
    for k = (N-1):-1:1
        x(k) = (b(k) - R(k,(k+1):N)*x((k+1):N))/R(k,k);
    end
    
end