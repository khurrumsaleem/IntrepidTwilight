classdef Shape < handle
    
    properties( SetAccess = protected , Hidden )
        Type = 'Shape'
    end
    
    properties
        
        % 2D location of shape (relative to some origin)
        %   Expected to be an array of (x,y) pairs
        Top
        Bottom
        Left
        Right

        % 3D information of shape
        SurfaceArea
        Volume
        
        % Kinematic geometry information
        FlowArea
        FlowDirection

    end
    
    methods

        function S = Shape(n)
            if (nargin >= 1) && not(isempty(n))
                S(n) = Shape();
            end
        end
        
        
        function [] = SetVertices(S,V1,V2,V3,V4)
            V = [V1(:,1),V2(:,1),V3(:,1),V4(:,1)];
            S.Left   = min(V,[],2);
            S.Right  = max(V,[],2);
            
            V = [V1(:,2),V2(:,2),V3(:,2),V4(:,2)];
            S.Top    = max(V,[],2);
            S.Bottom = min(V,[],2);
        end

    end
    
end