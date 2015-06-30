function ie = ImplicitEuler(spatialDiscretization)

    %   The SimpleTemporalDiscretization() is implicit euler, but the closure exists
    %   to act as a base for other linear-multistep schemes.
    ie = IntrepidTwilight.TransientStride.SimpleTemporalDiscretization(spatialDiscretization);

end