classdef Shape < handle
    
    properties( SetAccess = protected , Hidden )
        Type = 'Shape'
        ID   = 'Intrepid:AdamantWave:Shape';
    end
    
    properties
        
        % 2D location of shape vertices (relative to some origin)
        %   The n-th vertex of the m-th shape lies in the n-th column
        %   of the m-th row.
        Vx
        Vy
        
        % 3D information of shape
        SurfaceArea
        Volume
        
    end
    
    methods
        
        function S = Shape(n)
            if (nargin >= 1) && not(isempty(n))
                S(n) = Shape();
            end
        end
        
        
        function [] = SetVertices(S,varargin)
            
            % Vertices
            Vertices  = cell2mat(varargin);
            VerticesX = Vertices(:,1:2:end);
            VerticesY = Vertices(:,2:2:end);
            
            % Count
            [Nshapes,Nvertices] = size(VerticesX);
            
            % Allocate instance storage (same end-information but unsorted)
            S.Vx = VerticesX;
            S.Vy = VerticesY;
            
            % Plug-in vertices sorted in ccw order
            for k = 1:Nshapes
                
                kConvex = S.RobustConvexHull(VerticesX(k,:),VerticesY(k,:));
                
                if length(kConvex)-1 ~= Nvertices
                    error([S.ID,':NonconvexShape'],...
                        'The specified vertices do not form a convex shape');
                end
                
                S.Vx(k,:) = VerticesX(k,kConvex(1:end-1));
                S.Vy(k,:) = VerticesY(k,kConvex(1:end-1));
                
            end
            
            
        end
        
    end
    
    
    
    methods(Static)
        function kConvexNew = RobustConvexHull(Vx,Vy)
            
            if isrow(Vx)
                Vx = Vx(:);
                Vy = Vy(:);
                Columnified = true;
            else
                Columnified = false;
            end
            
            % Set up for loop
            Npoints = length(Vx);
            
            %	Indices for vertices known to be on hull
            kConvex = convhull(Vx,Vy)    ;
            Nconvex = length(kConvex) - 1;
            
            
            %   Check if any vertices were determined to not be on the
            %   hull.  If there are, perform a robust check to make
            %   sure there are no collinear false-positives to within
            %   twice machine-epsilon.
            if Nconvex ~= Npoints
                
                %   Define supposedly non-convex points
                kNotConvex   = setdiff(1:Npoints,kConvex(1:Nconvex));
                P            = [Vx(kNotConvex),Vy(kNotConvex)];
                
                %   Allocation for insertion information
                kConvexNew  = zeros(Npoints,1);
                Bottom      = 2;
                Top         = 1;
                
                % Sweep through convex hull
                kConvexNew(1) = kConvex(1);
                for k = 1:Nconvex
                    
                    P1     = [Vx(kConvex(k))  ,Vy(kConvex(k))  ];
                    P2     = [Vx(kConvex(k+1)),Vy(kConvex(k+1))];
                    OnHull = IsCollinear(P1,P2,P,[],true);
                    
                    if any(OnHull)
                        
                        %   Order Points on the hull based on distance from P1
                        DistFromP1 = (P(OnHull,1) - P1(1)).^2  +...
                            (P(OnHull,2) - P1(2)).^2  ;
                        [~,Ordered] = sort(DistFromP1);
                        
                        % Insert into new convex index vector
                        Top = Top + nnz(OnHull);
                        kOnHull = kNotConvex(OnHull);
                        kConvexNew(Bottom:Top) = kOnHull(Ordered);
                        Bottom = Top + 1;
                        
                        % Contract outstanding vertices
                        P(OnHull,:) = [];
                        kNotConvex(OnHull) = [];
                        
                        if isempty(kNotConvex)
                            Top = Nconvex - (k+1) + Bottom;
                            kConvexNew(Bottom:Top) = kConvex(k+1:Nconvex);
                            break;
                        end
                    end
                    
                    % Add second vertex to the list
                    kConvexNew(Bottom) = kConvex(k+1);
                    Top    = Top + 1;
                    Bottom = Bottom + 1;
                    
                end
                
                kConvexNew = [kConvexNew(1:Top);kConvexNew(1)];
                
                if Columnified
                    kConvexNew = kConvexNew';
                end
                
            end
        end
        
    end
    
end


function TrueFalse = IsCollinear(P1,P2,P,Tolerance,CheckBetween)
    
    if (nargin < 4) || isempty(Tolerance)
        Tolerance = 2*eps();
    end
    
    if (nargin < 5) || isempty(CheckBetween)
        CheckBetween = false;
    end
    
    if size(P,2) ~= 2
        
        P1 = P1';
        P2 = P2';
        P  = P' ;
        
    end
    
    % X points
    x1 = P1(:,1);
    x2 = P2(:,1);
    x  = P (:,1);
    
    % Y points
    y1 = P1(:,2);
    y2 = P2(:,2);
    y  = P (:,2);
    
    % Colinear error
    yLine = (y2-y1)./(x2-x1).*(x - x1) + y1;
    Error = abs(yLine - y);
    
    % Vertical Line check
    IsVertCollinear = (x1 == x2) & (x1 == x);
    
    % Output
    TrueFalse = (Error < Tolerance) | IsVertCollinear;
    
    if CheckBetween
        x12       = [x1,x2];
        IsBetween = (min(x12,[],2) <= x) & (x <= max(x12,[],2));
        TrueFalse = TrueFalse & IsBetween;
    end
    
end