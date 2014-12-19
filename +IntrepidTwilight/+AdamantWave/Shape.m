classdef Shape < IntrepidTwilight.AdamantWave.AdamantParent
    
    properties(SetAccess = private , Hidden)
        ID   = 'Shape';
        Purpose = ...
            'Holds all geometric data and location information of a shape.';
    end
    
    properties
        % 2D Location of shape (relative to some origin)
        Top
        Bottom
        Left
        Right

        % 3D information of shape
        SurfaceArea
        Volume
        
        % Description of shape
        Description
    end
    
    methods
        function S = Shape(n)
            if (nargin >= 1) && not(isempty(n))
                S(n) = FlowInterface();
            end
        end
    end
    
end