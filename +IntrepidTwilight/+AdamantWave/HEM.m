function hem = HEM()


    %   Attach evolver instance directly to simulation base
    components.model               = IntrepidTwilight.AdamantWave.BasicModel();
    components.spacediscretization = 0;
    components.timediscretization  = 0;
    components.residual            = 0;
    components.preconditioner      = 0;
    components.solver              = 0;
    components.evolver             = 0;
    
    hem.set          = @(component,varargin) set(component,varargin{:}) ;
    hem.bind         = @(component) bind(component)                     ;
    hem.get          = @(component) get(component)                      ;
    hem.modelValue   = @(varargin) set('model'  ,varargin{:})           ;
    hem.evolverValue = @(varargin) set('evolver',varargin{:})           ;
    


    %   User options
    hem.model                = IntrepidTwilight.AdamantWave.BasicModel();
    hem.discretization.space = 'Quasi2DUpwind'                          ;
    hem.discretization.time  = 'ImplicitEuler'                          ;
    hem.solver.name          = 'JFNK'                                   ;
    hem.solver.scheme        = 'two-level'                              ;
    hem.preconditioner       = 'none'                                   ;



    
    %   Finalizes options and builds components
    hem.build  = @() build()    ;
    hem.evolve = @() evolve()   ;
    hem.run    = @() evolve()   ;
    hem.data   = @() getData()  ;


    
    
    function [] = build()
        
        
        switch(lower(hem.solver.scheme))

            case('two-level')                
                hem = IntrepidTwilight.AdamantWave.HEM_TwoLevel(hem);



            case('coupled')
                args = {hem.discretization.space,components.model,struct('scheme','coupled')};
                components.spacediscretization = buildComponent('AdamantWave',args{:});
                %
                %   Build temporal discretization
                components.timediscretization = buildComponent('TransientStride',hem.discretization.time,components.spacediscretization);
                %
                %   Build residual
                components.residual = buildComponent('executive','Residual',components.timediscretization);
                %
                %   Build preconditioner
                components.preconditioner = buildComponent('executive','Preconditioner',components.residual,components.preconditioner);
                %
                %   Build solver
                components.solver = buildComponent('TenaciousReduction',hem.solver.name,components.residual,components.preconditioner);
                %
                %   Build evolver
                components.evolver = buildComponent('executive','Evolver',components.solver,components.residual);

            otherwise
        end
        
    end
    
    
    function [] = evolve()
        components.evolver.evolve();
    end
    
    function data = getData()
        data = components.evolver.getData();
    end
    function [] = bind(object)
        if isfield(components,object.type)
            components.(object.type) = object;
        end
    end
    function [] = set(component,varargin)
        components.(component).set(varargin{:});
    end
    function comp = get(component)
        comp = components.(component).get();
    end

end

function component = buildComponent(module,object,varargin)
    component = IntrepidTwilight.executive.build(module,object,varargin{:});
end
