classdef FlowInterface < IntrepidTwilight.AdamantWave.AdamantParent
    
    properties(SetAccess = private , Hidden)
        ID   = 'FlowInterface';
        Purpose = [...
            'Describes the necessary information to geometrically '...
            'define a fluid-fluid interface that permits flow.'];
    end
    
    properties
        Upwind
        Downwind
        Normal
        FlowArea
    end
    
    methods
        function FI = FlowInterface(n)
            if (nargin >= 1) && not(isempty(n))
                FI(n) = FlowInterface();
            end
        end
    end
    
end