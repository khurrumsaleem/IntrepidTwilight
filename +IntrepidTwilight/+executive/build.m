function component = build(module,objectName,varargin)

    objectRaw  = ['IntrepidTwilight.',module,'.',objectName];
    objectName = which(objectRaw)                           ;

    if not(isempty(objectName))

        objectName = regexpi(objectName,'\\([a-z0-9\_]+?)\.m','tokens');
        objectName = char(objectName{1});
        component  = IntrepidTwilight.(module).(objectName)(varargin{:});

    else
        
        error('IntrepidTwilight:executive:build:objectNotFound',...
            'The object ''%s'' could not be found.',objectRaw);

    end

end