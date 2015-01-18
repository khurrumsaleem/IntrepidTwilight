function Cinter = GetInterfaceIndices(iUp,iDown)
    
    % Order the inputs
    nInter = length(iUp)            ;
    nMC    = max([iUp(:);iDown(:)]) ;
    iInter = (1:length(iUp))'       ;
    one    = iInter./iInter         ;
    
    % Summation matrix
    Cinter = sparse(iUp(:)  ,iInter,-one,nMC,nInter) + ...
             sparse(iDown(:),iInter,+one,nMC,nInter);
    
    
end