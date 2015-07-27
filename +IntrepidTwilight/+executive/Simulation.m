function simulation = Simulation(moduleName)
    
    if (nargin >= 1) && not(isempty(moduleName))
        setBuildModuleName(moduleName);
    end
    
    
    %   Define the components of a Simulation
    componentNames = {'model','spacediscretization','timediscretization',...
        'residual','preconditioner','solver','evolver'};
    moduleNames    = {'','','TransientStride','executive','executive',...
        'TenaciousReduction','executive'};
    
    
    %   Allocate containers
    for component = componentNames 
        components.(component{1}) = []  ;   %   Container for object instances
        names.(component{1})      = []  ;   %   Container for object choices
        isNotBuilt.(component{1}) = true;   %   Build indicators
    end

    
    %   Public methods
    simulation.type   = 'simulation'                                    ;
    simulation.is     = @(s) strcmpi(s,simulation.type)                 ;
    simulation.bind   = @(component) bind(component)                    ;
    simulation.get    = @(component) get(component)                     ;
    simulation.build  = @(varargin) build(varargin{:})                  ;
    simulation.evolve = @() evolve()                                    ;
    simulation.run    = @() evolve()                                    ;
    simulation.data   = @() getData()                                   ;
    
    
    %   Generates a choice list for the passed instance
    simulation.makeChoice = @(instance,key,choices) makeChoice(instance,key,choices);

    %   Late module name bind
    simulation.setBuildModuleName = @(m) setBuildModuleName(m);
    
    
    %   Create specific component setters for convenience
    for component = componentNames
        simulation.(component{1}).set = @(key,value) components.(component{1}).set(key,value)   ;
        simulation.(component{1}).get = @()          components.(component{1})                  ;
    end
    
    
    
    
    % ========================================= %
    %            Methods Defintions             %
    % ========================================= %
    
    %   Evolve the solution
    function [] = evolve()
        
        if not(any(structfun(@(v)v,isNotBuilt)))
        
            %   Prepare components
            for componentName = componentNames
                if isfield(components.(componentName{1}),'prepare')
                    components.(componentName{1}).prepare();
                end
            end
            
            %   Execute evolution
            components.evolver.evolve();
            
        else
            
            error('IntrepidTwilight:executive:Simulation:NotAllBuilt',...
            'Simulation could not be evolved because not all components are built');
            
        end

    end
    
    %   Retrieve evolution data
    function data = getData()
        data = components.evolver.getData();
    end
    
    %   Bind object to components struct based on it declared type
    function [] = bind(object)
        if isstruct(object) && isfield(components,object.type)
            components.(object.type) = object   ;
            isNotBuilt.(object.type) = false    ;
        end
    end
    
    
    %   Set object data within given component
    function [] = setBuildModuleName(component,value)
        moduleNames.(component) = value;
    end
    
    
    
    
    
    function object = makeChoice(object,key,choices)
        
        choices = choices(:)';
        for choice = choices
            object.choose.(key).(choice{1}) = @() setComponentName(choice{1});
        end
        
        function [] = setComponentName(name)
            names.(key) = name;
        end
        
    end
    
    
    
    
    
    
    % =================================================== %
    %                    Build System                     %
    % =================================================== %
    
    function [] = build()
        
        
        if not(isempty(moduleName)) && not(all(structfun(@(n)isempty(n),names)))
            
            for component = componentNames  %#ok<*FXUP>
                if isNotBuilt.(component{1})
                    
                    switch(component{1})
                        case('model')
                            module    = moduleNames.model            ;
                            name      = names.spacediscretization   ;
                            arguments = {}                          ;
                            
                        case('spacediscretization')
                            module     = moduleNames.spacediscretization ;
                            name       = names.spacediscretization      ;
                            arguments  = {components.model}             ;
                            
                        case('timediscretization')
                            module     = 'TransientStride'                  ;
                            name       = names.timediscretization           ;
                            arguments  = {components.spacediscretization}   ;
                            
                        case('residual')
                            module     = 'executive'                    ;
                            name       = names.residual                 ;
                            arguments  = {components.timediscretization};
                            
                        case('preconditioner')
                            module     = 'executive'            ;
                            name       = names.preconditioner   ;
                            arguments  = {components.residual}  ;
                            
                        case('solver')
                            module     = 'TenaciousReduction'                           ;
                            name       = names.solver                                   ;
                            arguments  = {components.residual,components.preconditioner};
                            
                        case('evolver')
                            module     = 'executive'                            ;
                            name       = names.evolver                          ;
                            arguments  = {components.solver,components.residual};
                            
                        otherwise
                            error('IntrepidTwilight:executive:build:UnrecognizedComponent',...
                                'An unknown component ''%s'' build request.',...
                                component{1});
                    end
                    
                    
                    
                    components.(component{1}) = ...
                        IntrepidTwilight.executive.build(module,name,arguments{:});
                    
                    if isstruct(components.(component{1})) && components.(component{1}).is(component{1})
                        isNotBuilt.(component{1}) = false;
                    else
                        components.(component{1}) = [];
                    end
                    
                end
            end
            
        elseif isempty(moduleName)

            %   Module is not set
            error('IntrepidTwilight:executive:build:NoModuleNameDefined',...
                'Component ''%s'' could not be built because module is not specified',...
                component{1});

        else % not(all(structfun(@(n)isempty(n),names)))
            
            %   All components names must be given
            error('IntrepidTwilight:executive:build:NotAllObjectChosen',...
                'Component ''%s'' could not be built because not all object names are set',...
                component{1});
        end
        
        
    end
    
end