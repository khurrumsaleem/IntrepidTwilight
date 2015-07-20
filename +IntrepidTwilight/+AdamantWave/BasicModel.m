function model = BasicModel()
    
    
    values.dimensionalizer.mass     = [];
    values.dimensionalizer.energy   = [];
    values.dimensionalizer.momentum = [];

    
    % ========================================================= %
    %                       State Fields                        %
    % ========================================================= %

    %   Initial 
    values.controlVolume.mass            = [];
    values.controlVolume.energy          = [];
    values.controlVolume.pressure        = [];
    values.controlVolume.temperature     = [];
    values.controlVolume.internalEnergy  = [];
    values.controlVolume.enthalpy        = [];
    values.controlVolume.entropy         = [];
    values.controlVolume.volume          = [];
    values.controlVolume.source.mass     = [];
    values.controlVolume.source.energy   = [];


    % ========================================================= %
    %                     Geometry Fields                       %
    % ========================================================= %

    %   Control volume connections
    values.momentumCell.from = [];
    values.momentumCell.to   = [];
    
    
    %   Flow direction
    values.momentumCell.momentum        = [];
    values.momentumCell.directionX      = [];
    values.momentumCell.directionY      = [];
    values.momentumCell.volumeFrom      = [];
    values.momentumCell.volumeTo        = [];
    values.momentumCell.LoD             = [];
    values.momentumCell.source.momentum = [];
    values.momentumCell.source.friction = [];


    %   Interfaces
    values.interface.up       = [];  %   Upwind volume
    values.interface.down     = [];  %   Downwind volume
    values.interface.normalX  = [];  %   Surface normal
    values.interface.normalY  = [];  %   Surface normal
    values.interface.flowArea = [];  %   Flow area
    
    
    

    model.type = 'model';
    model.is   = @(s) strcmpi(s,model.type);
    model.set  = @(varargin) set(varargin{:});
    model.get  = @(varargin) get(varargin{:});
    
    function [] = set(varargin)
        values = setfield(values,varargin{:});
    end
    function value = get(varargin)
        if (nargin == 0)
            value = values;
        else
            value = getfield(values,varargin{:});
        end
    end
end
