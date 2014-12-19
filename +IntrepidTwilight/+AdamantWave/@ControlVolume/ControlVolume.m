classdef ControlVolume < IntrepidTwilight.AdamantWave.AdamantParent
    
    properties(SetAccess = private , Hidden)
        ID   = 'ControlVolume';
        Purpose = ...
            'To describe the state of a (finite) fluid-filled volume.';
    end
    
    properties
        Geometry = IntrepidTwilight.AdamantWave.Shape();
        State    = IntrepidTwilight.AdamantWave.State();
    end
    
    
    methods
        function CV = ControlVolume(n)
            if (nargin >= 1) && not(isempty(n))
                CV(n) = IntrepidTwilight.AdamantWave.ControlVolume();
            end
        end
    end
    
end