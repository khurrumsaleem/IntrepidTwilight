classdef State < IntrepidTwilight.AdamantWave.AdamantParent
    
    properties(SetAccess = private , Hidden)
        ID   = 'State';
        Purpose = ...
            'Store kinematic and thermodynamic data.';
    end
    
    properties
        Fluid

        % Thermodynamic properties
        Pressure
        Temperature
        Density
        Volume
        Entropy
        Enthalpy

        % Transport properties
        ThermalConductivity
        DynamicViscosity

        % Kinemtica properties
        Momentum
        Velocity
        
    end
    
    methods
        function S = State(n)
            if (nargin >= 1) && not(isempty(n))
                S(n) = State();
            end
        end
    end
    
end