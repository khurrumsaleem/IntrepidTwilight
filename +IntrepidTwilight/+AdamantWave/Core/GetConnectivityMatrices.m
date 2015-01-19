function [Ccv,Cmc] = GetConnectivityMatrices(fromTo)
%
%   Build the connectivity matrices for both control volumes (Ccv) and
%   momentum cells (Cmc) given a list of "from" and "to" values to build
%   the associated graphs.
%
%   Currently, this function is NOT so-called ID "safe".  This means that
%   the sorted, unique entries of the FromTo input are expected to 
%   monotonically run from 1 to Ncv by 1s.  This will be fixed in the
%   future.
%

    % Semantically partition
    Froms = fromTo(:,1);
    Tos   = fromTo(:,2);


    % Determine the number of momentum cells
    Nmc = size(fromTo,1);
    
    % Dtermine the number of control volumes
    Uniques = unique(fromTo(:));
    Ncv     = length(Uniques);
    
    
    % ================================================================== %
    %                       Momentum Cell Matrix                         %
    % ================================================================== %
    %   Create this first since it is easier.
    I   = 1:Nmc;
    Cmc = sparse(Ncv,Nmc)                   ;
    Cmc = sparse(Froms,I,-1,Ncv,Nmc) + Cmc  ;
    Cmc = sparse(Tos  ,I,+1,Ncv,Nmc) + Cmc  ;
    

    
    % ================================================================== %
    %                       Control Volume Matrix                        %
    % ================================================================== %
    %   Get the number of unique control volumes for allocation
    
    %   Since control volumes are considered connected to themselves for 
    %   donoring, begin the build with a sparse identity matrix.
    Ccv = speye(Ncv);
    
    for k = 1:Ncv
        CV           = Uniques(k)                               ;
        ConnectedCVs = [Tos(Froms == CV);Froms(Tos == CV)]      ;
        Ccv          = sparse(CV,ConnectedCVs,1,Ncv,Ncv) + Ccv  ;
    end
 
end