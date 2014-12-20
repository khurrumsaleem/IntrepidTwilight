classdef Kinematics < handle

    properties
        
        % Speeds
        Momentum
        Velocity
        MassFlowRate
        
    end
    
    methods
        function Kin = Kinematics(n)
            if (nargin >= 1) && not(isempty(n))
                Kin(n) = Kinematics();
            end
        end
    end
    
end