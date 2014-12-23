classdef Shape < handle
    
    properties( SetAccess = protected , Hidden )
        Type = 'Shape'
    end
    
    properties
        
        % 2D location of shape (relative to some origin)
        %   Expected to be an array of (x,y,z,...) tuples
        V

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
            
            Nvert = length(varargin);
            V     = cell2mat(varargin);
            
            Vx   = V(:,1:2:end);
            Vy   = V(:,2:2:end);


            R    = Vx.^2 + Vy.^2;
            Rmin = min(R,[],2);
            
            S.V1
            
            
            
        end

    end
    
end