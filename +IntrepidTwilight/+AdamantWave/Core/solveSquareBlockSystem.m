function x = solveSquareBlockSystem(D,b)
    
    
    [nUnknowns,blockWidth] = size(D);
    nBlocks = nUnknowns/blockWidth;
    
    if (nBlocks ~= round(nBlocks))
        error('solveSquareBlockSystem:NonSquareBlockStructure',...
            ['Compressed block input ''D'' does not have an integer',...
            'number of blocks']);
    end
    
    x = b;
    for k = 1:nBlocks
        iBlock = ((k-1)*blockWidth+1):(k*blockWidth);
        x(iBlock) = D(iBlock,:)\b(iBlock);
    end
    
end