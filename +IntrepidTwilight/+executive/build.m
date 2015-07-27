function component = build(module,objectName,varargin)

    objectName = which(['IntrepidTwilight.',module,'.',objectName]);

    if not(isempty(objectName))
        objectName = regexpi(objectName,'\\([a-z0-9\_]+?)\.m','tokens');
        objectName = char(objectName{1});
        component  = IntrepidTwilight.(module).(objectName)(varargin{:});

    else
        component = [];

    end

end