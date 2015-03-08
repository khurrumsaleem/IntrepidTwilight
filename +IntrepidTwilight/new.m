function object = new(request,varargin)
    
    
    switch(lower(request))
        case('problem')
            object = IntrepidTwilight.executive.Problem();


        case('simulation')
            if (length(varargin) >= 1)
                object = IntrepidTwilight.executive.Simulation(varargin{1});
            else
                error('IntrepidTwilight:new:noProblemData',...
                    'Required input problem data no present for new simulation.');
            end


        otherwise
            error('IntrepidTwilight:new:unsupportedRequest',...
                'Submitted unknown request ''%s''.',request);
    end
    
    
end