function hem = HEM()


    %   Closure variables
    spacediscretization = 0;
    timediscretization  = 0;
    solver              = 0;
    residual            = 0;
    preconditioner      = 0;
    simulation          = 0;


    %   User Options
    hem.model                = IntrepidTwilight.AdamantWave.BasicModel();
    hem.discretization.space = 'Quasi2DUpwind'                          ;
    hem.discretization.time  = 'ImplicitEuler'                          ;
    hem.solver.name          = 'JFNK'                                   ;
    hem.solver.scheme        = 'segregated'                             ;
    hem.preconditioner       = 'none'                                   ;



    
    %   Finalizes options and builds components
    hem.processOptions = @() processOptions();


    
    
    function [] = processOptions()
        
        %   Build spatial discretization
        spacediscretization = buildComponent('AdamantWave',hem.discretization.space,hem.model);
        
        
        switch(lower(hem.solver.scheme))
            case('segregated')
            case('coupled')
            otherwise
        end



        %   Build temporal discretization
        timediscretization = buildComponent('TransientStride',hem.discretization.time,spacediscretization);
        
        %   Build residual
        residual = buildComponent('executive','Residual',timediscretization);
        
        %   Build preconditioner
        preconditioner = buildComponent('executive','Preconditioner',residual,hem.preconditioner);
        
        %   Build solver
        solver = buildComponent('TenaciousReduction',hem.solver.name,residual,preconditioner);
        
        %   Build simulation
        simulation = buildComponent('executive','Simulation',solver,residual);
        
        
    end

end

function component = buildComponent(module,object,varargin)
    objectName = which(['IntrepidTwilight.',module,'.',object]);
    if not(isempty(objectName))
        objectName = regexpi(objectName,'\\([a-z0-9\_]+?)\.m','tokens');
        objectName = char(objectName{1});
        component  = IntrepidTwilight.(module).(objectName)(varargin{:});
    else
        component = [];
    end
end