function sd = SpaceDiscretization(model)

    %   Inherit and overwrite ID information
    sd         = IntrepidTwilight.executive.Object()   ;
    sd.type    = 'spaceDiscretization'                 ;
    sd.is      = @(s) strcmpi(s,sd.type)               ;
    
    
    %   Bind/Retrieve model value
    sd.bind     = @(m) bind(m)  ;
    sd.retrieve = @() retrieve();
    %
    function [] = bind(m)
        model = m;
    end
    function m = retrieve()
        m = model;
    end



    %   Handle argument
    if (nargin > 0) && isstruct(model) && model.is('model');
        bind(model);
    else
        bind([]);
    end

end