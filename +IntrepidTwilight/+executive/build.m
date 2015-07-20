function component = build(module,object,varargin)

    objectName = which(['IntrepidTwilight.',module,'.',object]);

    if not(isempty(objectName))
        objectName = regexpi(objectName,'\\([a-z0-9\_]+?)\.m','tokens');
        objectName = char(objectName{1});
        component  = IntrepidTwilight.(module).(objectName)(varargin{:});

    else
        component = [];

    end

end