function x = GMRES(Afun,b,x0,preRight)
    
    
    %   Runtime parameters
    nu        = 0.90        ;
    tolerance = 1E-10       ;
    n         = numel(x0)   ;
    m         = n           ;
    r0        = b - Afun(x0);
    r0Norm    = norm(r0,2)  ;
    
    %   Allocation
    Zeros = zeros(n,1)          ;
    Z     = zeros(n,m)          ;
    R     = Z                   ;
    H     = Z                   ;
    e     = [1 ; Zeros(1:n-1)]  ;
    alpha = e                   ;
    
    
    %   First basis vector for update and A*z1 product
    Z(:,1) = r0 / r0Norm            ;
    R(:,1) = Afun(preRight(Z(:,1))) ;
    
    % Compute Householder vector and bring R(:,1) into upper triangular form
    h      = R(:,1)                                 ;
    h      = -Signum(h(1)) * norm(h,2) * e - h      ;
    H(:,1) = h / norm(h,2)                          ;
    R(:,1) = R(:,1) - 2 * H(:,1) * (H(:,1)'*R(:,1)) ;
    
    
    % Get the first column of the unitary matrix
    q = e - 2 * H(:,1) * (H(:,1)'*e);
    e = e(1:n-1)                    ;
    
    % Residual update
    alpha(1) = q'*r0            ;
    rk       = r0 - alpha(1)*q  ;
    
    % Assign residual norms to determine which basis to use
    k        = 1            ;
    rkm1Norm = r0Norm       ;
    rkNorm   = norm(rk,2)   ;
    
    if rkNorm > tolerance
        for k = 2:m
            
            % Choose the next basis vector
            if rkNorm <= nu*rkm1Norm
                Z(:,k) = rk/rkNorm  ;   %   GCR (RB-SGMRES) basis
            else
                Z(:,k) = q          ;   %   Simpler GMRES basis
            end
            
            % Compute, store A*zk in R and apply all previous projections to new the column
            R(:,k) = Afun(preRight(Z(:,k)));
            for m = 1:k-1
                R(:,k) = R(:,k) - H(:,m)*(2*H(:,m)'*R(:,k));
            end
            
            % Get the next Householder vector
            h        = R(k:n,k)                     ;
            h        = -Signum(h(1))*norm(h,2)*e - h;
            h        = h ./ norm(h,2)               ;
            H(k:n,k) = h                            ;
            
            %   Apply new projection to R to bring it into upper triangular form;
            for m = 1:k
                R(:,m) = R(:,m) - 2 * H(:,k) * (H(:,k)'*R(:,m));
            end
            
            
            % Get the k-th column of the current unitary matrix
            q = [Zeros(1:k-1) ; e - 2*h*(h'*e)];
            for m = k-1:-1:1
                q = q - 2*H(:,m)*(H(:,m)'*q);
            end
            e = e(1:end-1);
            
            % Update residual
            alpha(k) = q'*rk            ;
            rk       = rk - alpha(k)*q  ;
            
            % Update residual norms
            rkm1Norm = rkNorm;
            rkNorm   = norm(rk,2);
            
            % Solve least-squares problem
            if rkNorm < tolerance
                break;
            end
            
        end
    end
    
    % Update to x
    Rtilde = triu(R(1:k,1:k));
    
    if rcond(Rtilde) > 100*eps()
        yk = Rtilde \ alpha(1:k);   % Solve the least-squares problem
    else
        %   Attempt to re-scale
        S  = diag(1./diag(Rtilde))  ;
        yk = (Rtilde*S) \ alpha(1:k);
        yk = S*yk                   ;
    end
    
    x = Z(:,1:k) * yk   ; % Calculate solution vector
    x = preRight(x)     ;
    
    
end

function s = Signum(s)
    if (s ~= 0)
        s = sign(s);
    else
        s = 1;
    end
end

