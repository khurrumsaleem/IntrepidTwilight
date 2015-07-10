function model = BasicModel()

    
    % ========================================================= %
    %                       State Fields                        %
    % ========================================================= %

    %   Initial 
    model.controlVolume.mass            = [];
    model.controlVolume.energy          = [];
    model.controlVolume.pressure        = [];
    model.controlVolume.temperature     = [];
    model.controlVolume.internalEnergy  = [];
    model.controlVolume.enthalpy        = [];
    model.controlVolume.entropy         = [];
    model.controlVolume.volume          = [];
    model.controlVolume.source.mass     = [];
    model.controlVolume.source.energy   = [];


    % ========================================================= %
    %                     Geometry Fields                       %
    % ========================================================= %

    %   Control volume connections
    model.momentumCell.from = [];
    model.momentumCell.to   = [];
    
    
    %   Flow direction
    model.momentumCell.momentum        = [];
    model.momentumCell.directionX      = [];
    model.momentumCell.directionY      = [];
    model.momentumCell.volumeFrom      = [];
    model.momentumCell.volumeTo        = [];
    model.momentumCell.LoD             = [];
    model.momentumCell.source.momentum = [];


    %   Interfaces
    model.interface.up       = [];  %   Upwind volume
    model.interface.down     = [];  %   Downwind volume
    model.interface.normalX  = [];  %   Surface normal
    model.interface.normalY  = [];  %   Surface normal
    model.interface.flowArea = [];  %   Flow area

end