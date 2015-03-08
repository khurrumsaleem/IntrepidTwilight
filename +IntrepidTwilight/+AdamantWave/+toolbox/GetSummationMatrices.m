function [Ccv,Cmc,Cp,inter] = GetSummationMatrices(fromTo,upDown,zDotN)

    % Semantically partition
    from     = fromTo(:,1)  ;
    to       = fromTo(:,2)  ;
    up       = upDown(:,1)  ;
    down     = upDown(:,2)  ;
    upDotN   = zDotN(:,1)   ;
    downDotN = zDotN(:,2)   ;

    % Determine the number of momentum cells
    nMC = size(fromTo,1);
    
    % Dtermine the number of control volumes
    Uniques = unique(fromTo(:));
    nCV     = length(Uniques);
    
    
    %   Sums over all MCs connected to CVs
    I   = 1:nMC;
    Ccv = sparse(nCV,nMC)                  ;
    Ccv = sparse(from,I,-1,nCV,nMC) + Ccv  ;
    Ccv = sparse(to  ,I,+1,nCV,nMC) + Ccv  ;



    %   Order the inputs
    nInter = length(up) ;
    iInter = (1:nInter)';

    % Sums over all interfaces between two MCs for the rho v^2 term
    Cmc = sparse(up  ,iInter,-1,nMC,nInter) + ...
          sparse(down,iInter,+1,nMC,nInter) ;
    
      
      
      
    % Sums over all interfaces between two MCs  for the P term
    iInter = zeros(2*nInter,1)  ;
    jInter = zeros(2*nInter,1)  ;
    sInter = zeros(2*nInter,1)  ;
    inter  = zeros(nInter,1)    ;
    for k = 1:nInter
        % Up
        iInter(k) = up(k)    ;
        jInter(k) = k        ;
        sInter(k) = upDotN(k);

        % Down
        iInter(k+nInter) = down(k)      ;
        jInter(k+nInter) = k            ;
        sInter(k+nInter) = -downDotN(k) ;
        
        % Interface control volume
            inter(k) = intersect(fromTo(up(k),:),fromTo(down(k),:));
    end
    Cp = sparse(iInter,jInter,sInter);
    
 
end