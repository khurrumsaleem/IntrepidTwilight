classdef ControlVolume < handle
    
    properties
        
        Name
        ID
        
        Geometry       = IntrepidTwilight.AdamantWave.Shape()           ;
        Thermodynamics = IntrepidTwilight.AdamantWave.Thermodynamics()  ;
        Kinematics     = IntrepidTwilight.AdamantWave.Kinematics()      ;
    end
    
    methods
        function CV = ControlVolume(n)
            if (nargin >= 1) && not(isempty(n))
                CV(n) = IntrepidTwilight.AdamantWave.ControlVolume();
            end
        end
    end
    
end