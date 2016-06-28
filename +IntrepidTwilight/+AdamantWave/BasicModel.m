function model = BasicModel(config)
    
    model = IntrepidTwilight.executive.Component();
    model = model.changeID(model,'BasicModel','model');
    
    
    model.set('dimensionalizer.mass'    , []);
    model.set('dimensionalizer.energy'  , []);
    model.set('dimensionalizer.momentum', []);
    model.set('dimensionalizer.volume'  , []);

    
    % ========================================================= %
    %                       State Fields                        %
    % ========================================================= %

    %   Initial 
    model.set('controlVolume.mass'          , []    );
    model.set('controlVolume.energy'        , []    );
    model.set('controlVolume.pressure'      , []    );
    model.set('controlVolume.temperature'   , []    );
    model.set('controlVolume.internalEnergy', []    );
    model.set('controlVolume.enthalpy'      , []    );
    model.set('controlVolume.entropy'       , []    );
    model.set('controlVolume.volume'        , []    );
    model.set('controlVolume.isElastic'     , []    );
    model.set('controlVolume.source.mass'   , @()[] );
    model.set('controlVolume.source.energy' , @()[] );


    % ========================================================= %
    %                     Geometry Fields                       %
    % ========================================================= %

    %   Control volume connections
    model.set('momentumCell.from', []);
    model.set('momentumCell.to'  , []);
    
    
    %   Flow direction
    model.set('momentumCell.momentum'           , []    );
    model.set('momentumCell.directionX'         , []    );
    model.set('momentumCell.directionY'         , []    );
    model.set('momentumCell.volumeFractionFrom' , []    );
    model.set('momentumCell.volumeFractionTo'   , []    );
    model.set('momentumCell.LoD'                , []    );
    model.set('momentumCell.source.momentum'    , @()[] );
    model.set('momentumCell.source.friction'    , []    );


    %   Interfaces
    model.set('interface.up'      , []);  %   Upwind volume
    model.set('interface.down'    , []);  %   Downwind volume
    model.set('interface.normalX' , []);  %   Surface normal
    model.set('interface.normalY' , []);  %   Surface normal
    model.set('interface.flowArea', []);  %   Flow area
    
    if (nargin > 0) && isstruct(config)
        model.set(config);
    end

end
