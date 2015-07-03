function model = BasicModel()

    
    % ========================================================= %
    %                       State Fields                        %
    % ========================================================= %

    %   Initial 
    model.controlVolume.mass           = [];
    model.controlVolume.energy         = [];
    model.controlVolume.pressure       = [];
    model.controlVolume.temperature    = [];
    model.controlVolume.internalEnergy = [];
    model.controlVolume.enthalpy       = [];
    model.controlVolume.entropy        = [];
    model.controlVolume.source.mass    = [];
    model.controlVolume.source.energy  = [];

    
    % ========================================================= %
    %                     Geometry Fields                       %
    % ========================================================= %

    %   Control volume connections
    model.momentumCell.from = [];
    model.momentumCell.to   = [];
    
    
    %   Flow direction
    model.momentumCell.directionX      = [];
    model.momentumCell.directionY      = [];
    model.momentumCell.volumeBack      = [];
    model.momentumCell.volumeFront     = [];
    model.momentumCell.LoD             = [];
    model.momentumCell.loss.friction   = [];
    model.momentumCell.source.momentum = [];

    
    %   Interfaces
    model.interface.volumeUp    = [];  %   Upwind volume
    model.interface.volumeDown  = [];  %   Downwind volume
    model.interface.normalX     = [];  %   Surface normal
    model.interface.normalY     = [];  %   Surface normal
    model.interface.surfaceArea = [];  %   Flow area
    model.interface.loss.form   = [];  %   Form loss coefficient


end