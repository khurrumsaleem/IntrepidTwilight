function J = numericalJacobianBlockDiagonal(fun,x,xBlockSize,epsilon)
    
    narginchk(3,4);
    
    % ============================================================= %
    %                             Check f                           %
    % ============================================================= %
    extract = @(x,indices) subsref(x,struct('type','()','subs',{{indices}}));
    
    if isa(fun,'function_handle')

        if (nargin(fun) >= 2)
            f    = @(x,k,ind) fun(x,k)  ;
            fOld = f                    ;
        elseif (nargin(fun) == 1)
            f    = @(x,k,ind) extract(fun(x),ind)   ;
            f0   = fun(x)                           ;
            fOld = @(x,k,ind) f0(ind)               ;
        else
            error('Function handle ''f'' must except at least 1 argument.');
        end

    elseif iscell(fun)

        if isCellFuntion(fun)

            if (nargin(fun{1}) >= 1)
                f    = @(x,k,ind) fun{k}(x) ;
                fOld = f                    ; 
            else
                error('Function handle ''f'' must except at least 1 argument.');
            end

        else
            error('Cell array ''f'' must be a collection of function handles.');
        end

    else
        error('Input ''f'' should be a function handle or cell array of function handles.');
    end




    % ============================================================= %
    %                            Check x0                           %
    % ============================================================= %
    if not(isnumeric(x))
        error('Input ''x0'' must be a numeric vector');
    end


    % ============================================================= %
    %                        Check xblockSize                       %
    % ============================================================= %
    nX = numel(x);
    if isscalar(xBlockSize)
        columns  = xBlockSize                   ;
        kBlocks  = nX/xBlockSize                ;
        rows     = xBlockSize * kBlocks         ;
        mColumns = xBlockSize(ones(1,kBlocks))  ;
    else
        columns  = max(xBlockSize)  ;
        kBlocks  = numel(xBlockSize);
        rows     = sum(xBlockSize)  ;
        mColumns = xBlockSize       ;
    end

    
    J(rows,columns) = 0         ;
    e               = (1:rows)' ;
    
    bottom = 1              ;
    top    = mColumns(1)    ;
    I      = bottom:top     ;
    for k = 1:kBlocks
        
        f0 = fOld(x,k,I);

        for m = 1:mColumns(k)

            J(I,m) = (f( x + (e == e(bottom + m - 1))*epsilon,k,I) - f0)/epsilon;
            
        end
        
        bottom = top + 1                            ;
        top    = top + mColumns( min(kBlocks,k+1) ) ;
        I      = bottom:top                         ;

    end
    
    
    
end

function TF = isCellFuntion(C)
    TF = all(cellfun(@(c) isa(c,'function_handle'),C));
end