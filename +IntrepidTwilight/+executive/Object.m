function object = Object()

    object.type = 'object'                      ;
    object.is   = @(s) strcmpi(s,object.type)   ;
    object.set  = @(key,value) set(key,value)   ;
    object.get  = @(key)       get(key)         ;
    
    parameters = struct();
    
    function [] = set(key,value)
        keys       = strsplit(key,'.')                  ;
        parameters = setfield(parameters,keys{:},value) ;
    end

    function value = get(key)
        keys  = strsplit(key,'.')           ;
        value = getfield(parameters,keys{:});
    end
    
end