function kappa = conditionNumberMaximum(A,p)

    narginchk(2,2);

    switch(class(A))

        case('double')
            kappa = cond(A,p);

        case('cell')
            kappa = max(cellfun(@(c) cond(c,p),A));

        otherwise
            error('IntrepidTwilight:ConvenientMeans:conditionNumberMaximum:invalidInputType',...
                'First input should be a double array or cell array of double arrays.'); 
    end

end