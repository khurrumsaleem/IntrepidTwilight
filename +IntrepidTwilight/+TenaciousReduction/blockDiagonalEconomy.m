function x = blockDiagonalEconomy(D,x,blockSize)
    
    x = mat2cell(x,blockSize);

    for k = 1:numel(D)
        x{k} = D{k} \ x{k};
    end

    x = vertcat(x{:});

end