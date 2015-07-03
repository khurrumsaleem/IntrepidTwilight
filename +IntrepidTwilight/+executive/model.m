function md = model()

    
    % ========================================================= %
    %                       State Fields                        %
    % ========================================================= %

    %   Initial 
    md.controlVolume.mass           = [];
    md.controlVolume.energy         = [];
    md.controlVolume.pressure       = [];
    md.controlVolume.temperature    = [];
    md.controlVolume.internalEnergy = [];
    md.controlVolume.enthalpy       = [];
    md.controlVolume.entropy        = [];
    md.controlVolume.source.mass    = [];
    md.controlVolume.source.energy  = [];

    
    % ========================================================= %
    %                     Geometry Fields                       %
    % ========================================================= %

    %   Control volume connections
    md.momentumCell.from = [];
    md.momentumCell.to   = [];
    
    
    %   Flow direction
    md.momentumCell.directionX      = [];
    md.momentumCell.directionY      = [];
    md.momentumCell.volumeBack      = [];
    md.momentumCell.volumeFront     = [];
    md.momentumCell.LoD             = [];
    md.momentumCell.loss.friction   = [];
    md.momentumCell.source.momentum = [];

    
    %   Interfaces
    md.interface.volumeUp    = [];  %   Upwind volume
    md.interface.volumeDown  = [];  %   Downwind volume
    md.interface.normalX     = [];  %   Surface normal
    md.interface.normalY     = [];  %   Surface normal
    md.interface.surfaceArea = [];  %   Flow area
    md.interface.loss.form   = [];  %   Form loss coefficient


end