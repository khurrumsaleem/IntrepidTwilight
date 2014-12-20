classdef Thermodynamics < handle

    properties
        
        % Fluid name
        Fluid
        Phase

        % Thermodynamic properties
        Pressure
        Temperature
        InternalEnergy
        Density
        SpecificVolume
        Entropy
        Enthalpy

        % Transport properties
        ThermalConductivity
        DynamicViscosity

    end
    
    methods
        function Thermo = Thermodynamics(n)
            if (nargin >= 1) && not(isempty(n))
                Thermo(n) = Thermodynamics();
            end
        end
    end
    
end