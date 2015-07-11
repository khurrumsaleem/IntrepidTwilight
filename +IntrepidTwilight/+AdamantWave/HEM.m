function hem = HEM()


    %   Attach evolver instance directly to simulation base
    components.model               = IntrepidTwilight.AdamantWave.BasicModel();
    components.spacediscretization = 0;
    components.timediscretization  = 0;
    components.residual            = 0;
    components.preconditioner      = 0;
    components.solver              = 0;
    components.evolver             = 0;
    
    hem.set          = @(component,varargin) set(component,varargin{:});
    hem.modelValue   = @(varargin) set('model'  ,varargin{:});
    hem.evolverValue = @(varargin) set('evolver',varargin{:});
    


    %   User options
    hem.model                = IntrepidTwilight.AdamantWave.BasicModel();
    hem.discretization.space = 'Quasi2DUpwind'                          ;
    hem.discretization.time  = 'ImplicitEuler'                          ;
    hem.solver.name          = 'JFNK'                                   ;
    hem.solver.scheme        = 'segregated'                             ;
    hem.preconditioner       = 'none'                                   ;



    
    %   Finalizes options and builds components
    hem.build  = @() build()    ;
    hem.evolve = @() evolve()   ;
    hem.run    = @() evolve()   ;
    hem.data   = @() getData()  ;


    
    
    function [] = build()
        
        
        switch(lower(hem.solver.scheme))
            case('segregated')
                
                %   Build full sd
                args = {hem.discretization.space,components.model.get(),struct('scheme','segregated')};
                components.spacediscretization = buildComponent('AdamantWave',args{:});
                %
                %   Build CV block
                sdCV.sd     = components.spacediscretization;
                sdCV.update = @(time) sdCV.sd.update(time)      ;
                sdCV.rhs    = @(qCV)  sdCV.sd.rhsMassEnergy(qCV);
                %
                %   Build MC block
                sdMC.sd     = components.spacediscretization;
                sdMC.update = @(time) sdMC.sd.update(time);
                sdMC.rhs    = @(qMC) sdMC.sd.rhsMomentum(qMC);
                %
                %   Build time steppers
                components.timediscretization    = buildComponent('TransientStride',hem.discretization.time,sdMC);
                components.timediscretization(2) = buildComponent('TransientStride',hem.discretization.time,sdCV);
                %
                %   Build residuals
                components.residual    = buildComponent('executive','Residual',components.timediscretization(1));
                components.residual(2) = buildComponent('executive','Residual',components.timediscretization(2));
                %
                %   Build preconditioners
                components.preconditioner    = buildComponent('executive','Preconditioner',components.residual(1),'none');
                components.preconditioner(2) = buildComponent('executive','Preconditioner',components.residual(2),'none');
                %
                %   Build solver
                components.solver = buildComponent('TenaciousReduction',hem.solver.name,components.residual,components.preconditioner);
                %
                %   Build evolver
                evolverResidual.update = @(value,time,step) arrayfun(@(s,e) s.update(value{e},time,step),components.residual,1:2);
                evolverResidual.is     = @(s) strcmpi(s,'residual');
                components.evolver = buildComponent('executive','Evolver',components.solver,evolverResidual);
                
                
            case('coupled')
                args = {hem.discretization.space,hem.model,struct('scheme','coupled')};
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
    
    function [] = set(component,varargin)
        components.(component).set(varargin{:});
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