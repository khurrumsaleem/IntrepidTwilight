function list = moduleRegistry(whichList)
    
    if (nargin < 1)
        whichList = 'all';
    end
    
    %   Registry Lists
    physicsList = {...
        'AdamantWave'       ;...
        'ConcordantHelm'    ;...
        'ResonantConduct'   };
    
    numericsList = {...
        'TenaciousReduction'    ;...
        'TransientStride'       };
    
    
    %   Return requested list
    switch(lower(whichList))
        case('all')
            list = [physicsList;numericsList];
            
        case('physics')
            list = physicsList;
            
        case('numerics')
            list = numericsList;
            
    end
    
    
end