function Cinter = GetInterfaceIndices(up,down)
    
    % Order the inputs
    nInter = length(up)            ;
    nMC    = max([up(:);down(:)]) ;
    iInter = (1:length(up))'       ;
    
    % Summation matrix
    Cinter = sparse(up(:)  ,iInter,-1,nMC,nInter) + ...
             sparse(down(:),iInter,+1,nMC,nInter);
    
    
end