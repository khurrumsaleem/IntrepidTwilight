function object = Object()

    %   ID system
    object.type     = 'object'                          ;
    object.name     = object.type                       ;
    object.is       = @(s) strcmpi(s,object.type)       ;
    object.named    = @(s) strcmpi(s,object.name)       ;
    object.changeID = @(varargin) changeID(varargin{:}) ;

    function instance = changeID(instance,name,type)
        instance.name  = name                           ;
        instance.named = @(s) strcmpi(s,instance.name)  ;
        
        if (nargin == 3)
            instance.type  = type                           ;
            instance.is    = @(s) strcmpi(s,instance.type)  ;
        end
    end
    
    
    
    

    %   Key-store system
    object.set  = @(key,varargin) set(key,varargin{:})  ;
    object.get  = @(varargin)  get(varargin{:})         ;
    
    store = struct();
    
    function [] = set(key,value)
        
        narginchk(1,2);
        
        if (nargin == 2)
            keys  = strsplit(key,'.')               ;
            store = setfield(store,keys{:},value)   ;
        else
            if isstruct(key)
                store = IntrepidTwilight.ConvenientMeans.structMerge(store,key);
            else
                error('IntrepidTwilight:executive:Object:tooFewInputs',...
                    'The one input form a ''set'' may only be used with input type ''struct''');
            end
        end

    end
    function value = get(key)
        if (nargin >= 1)
            keys  = strsplit(key,'.')       ;
            value = getfield(store,keys{:}) ;
        else
            value = store;
        end
    end
    
end