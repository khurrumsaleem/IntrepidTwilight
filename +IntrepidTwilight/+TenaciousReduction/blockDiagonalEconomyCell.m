function x = blockDiagonalEconomyCell(D,x)

    for k = 1:numel(D)
        x{k} = D{k} \ x{k};
    end

    x = vertcat(x{:});

end