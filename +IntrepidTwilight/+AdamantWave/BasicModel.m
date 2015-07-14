function model = BasicModel()
    
    
    model.dimensionalizer.mass     = [];
    model.dimensionalizer.energy   = [];
    model.dimensionalizer.momentum = [];

    
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
    model.momentumCell.source.friction = [];


    %   Interfaces
    model.interface.up       = [];  %   Upwind volume
    model.interface.down     = [];  %   Downwind volume
    model.interface.normalX  = [];  %   Surface normal
    model.interface.normalY  = [];  %   Surface normal
    model.interface.flowArea = [];  %   Flow area
    
    
    

    model.type = 'model';
    model.is   = @(s) strcmpi(s,model.type);
    model.set  = @(varargin) set(varargin{:});
    model.get  = @(varargin) get(varargin{:});
    
    function [] = set(varargin)
        model = setfield(model,varargin{:});
    end
    function value = get(varargin)
        if (nargin == 0)
            value = model;
        else
            value = getfield(model,varargin{:});
        end
    end
end
