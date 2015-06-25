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

    
    % ========================================================= %
    %                     Geometry Fields                       %
    % ========================================================= %

    %   Control volume connections
    md.momentumCell.from = [];
    md.momentumCell.to   = [];
    
    
    %   Flow direction
    md.momentumCell.flowX       = [];
    md.momentumCell.flowY       = [];
    md.momentumCell.volumeBack  = [];
    md.momentumCell.volumeFront = [];
    md.momentumCell.LoD         = [];

    
    %   Interfaces
    md.interface.volumeUp   = [];  %   Upwind volume
    md.interface.volumeDown = [];  %   Downwind volume
    md.interface.normalX    = [];  %   Surface normal
    md.interface.normalY    = [];  %   Surface normal
    md.interface.area       = [];  %   Flow area
    
    
    
    %   Private meta-information
    meta.type    = 'model'      ;
    meta.package = 'AdamantWave';
    
    md.meta = @(s) meta.(s);
    
end