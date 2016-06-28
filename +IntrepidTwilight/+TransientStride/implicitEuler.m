function ie = ImplicitEuler(varargin)

    %   The TimeDiscretization() is implicit euler, but the closure exists
    %   to act as a base for other linear-multistep schemes.
    ie = IntrepidTwilight.TransientStride.TimeDiscretization(varargin{:});
    ie = ie.changeID(ie,'impliciteuler','timediscretization');

end