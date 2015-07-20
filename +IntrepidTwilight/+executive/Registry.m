function registry = Registry()
    
    items     = struct();
    registry.items = @() get();
    
    registry.set       = @(k,v)      set(k,v)           ;
    registry.get       = @(varargin) get(varargin{:})   ;
    registry.overwrite = @(s)        overwrite(s)       ;
    registry.clear     = @()         clear()            ;
    
    
    function [] = set(key,value)
        keys  = strsplit(key,'.');
        items = setfield(items,{1},keys{:},value);
    end
    
    function value = get(key)
        if (nargin > 0)
            keys  = strsplit(key,'.');
            value = getfield(items,{1},keys{:});
        else
            value = items;
        end
    end

    function [] = overwrite(s)
        if isstruct(s)
            items = s;
        end
    end
    
    function [] = clear()
        items = [];
    end
    
    
end