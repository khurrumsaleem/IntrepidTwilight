function simulation = Simulation(moduleName)
    
    %   Inherit
    simulation = IntrepidTwilight.executive.Object();
    
    
    if (nargin >= 1) && not(isempty(moduleName))
        setBuildModuleName(moduleName);
    end
    
    
    %   Define the components of a Simulation
    componentNames = {'model','spacediscretization','timediscretization',...
        'residual','preconditioner','solver','evolver'};
    moduleNameDefaults = {'','','TransientStride','executive','executive',...
        'TenaciousReduction','executive'};
    
    
    %   Allocate containers
    for m = 1:numel(componentNames)  %#ok<*FXUP>
        componentName               = componentNames{m}     ;
        moduleNames.(componentName) = moduleNameDefaults{m} ;
        components.(componentName)  = []                    ;   %   Container for object instances
        choices.(componentName)     = []                    ;   %   Container for object choices
        isNotBuilt.(componentName)  = true                  ;   %   Build indicators
    end

    
    %   Public methods
    simulation.type   = 'simulation'                    ;
    simulation.name   = simulation.type                 ;
    simulation.is     = @(s) strcmpi(s,simulation.type) ;
    simulation.named  = @(s) strcmpi(s,simulation.name) ;
    simulation.bind   = @(component) bind(component)    ;
    simulation.get    = @(component) get(component)     ;
    simulation.build  = @(varargin) build(varargin{:})  ;
    simulation.evolve = @() evolve()                    ;
    simulation.run    = @() evolve()                    ;
    simulation.data   = @() getData()                   ;
    
    
    %   Generates a choice list for the passed instance
    simulation.makeChoice = @(instance,key,choices) makeChoice(instance,key,choices);

    %   Late module name bind
    simulation.setBuildModuleName = @(component,name) setBuildModuleName(component,name);
    
    
    %   Create specific component setters for convenience
    for componentName = componentNames
        simulation.(componentName{1}).set = @(key,value) set([componentName{1},'.',key],value)  ;
        simulation.(componentName{1}).get = @()          get(componentName{1})                  ;
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
    function [] = setBuildModuleName(component,name)
        moduleNames.(component) = name;
    end
    
    
    
    
    
    function object = makeChoice(object,key,choiceList)
        
        choiceList = choiceList(:)';
        for choice = choiceList
            object.choose.(key).(choice{1}) = @() setChoiceName(choice{1})  ;
            object.choice.(key)             = @() getChoice(key)            ;
        end
        
        function [] = setChoiceName(name)
            choices.(key) = name;
        end
        function choice = getChoice(key)
            if isfield(choices,key)
                choice = choices.(key);
            else
                choice = [];
            end
        end
        
    end
    
    
    
    
    
    
    % =================================================== %
    %                    Build System                     %
    % =================================================== %
    
    function [] = build(bindAtBuild)
        
        if (nargin < 0) || isempty(bindAtBuild)
            bindAtBuild = true;
        end
        
        
        if not(isempty(moduleNames)) && not(all(structfun(@(n)isempty(n),choices)))
            
            if bindAtBuild
                arguments = {...
                    {},...
                    {components.model},...
                    {components.spacediscretization},...
                    {components.timediscretization},...
                    {components.residual},...
                    {components.residual,components.preconditioner},...
                    {components.solver,components.residual}};
            else
                arguments = {{},{},{},{},{},{},{}};
            end
            
            
            
            for k = 1:numel(componentNames)
                
                componentName = componentNames{k};
                
                if isNotBuilt.(componentName)
                    module = moduleNames.(componentName);
                    name   = choices.(componentName)    ;
                    list   = arguments{k}               ;
                                      
                    component = ...
                        IntrepidTwilight.executive.build(module,name,list{:});
                    
                    if isstruct(componentName) && component.is(componentName)
                        components.(componentName) = component  ;
                        isNotBuilt.(componentName) = false      ;
                    else
                        components.(componentName) = [];
                    end
                    
                end
            end
            
        elseif isempty(moduleName)

            %   Module is not set
            error('IntrepidTwilight:executive:build:NoModuleNameDefined',...
                'Component ''%s'' could not be built because module is not specified',...
                componentName);

        else % not(all(structfun(@(n)isempty(n),names)))
            
            %   All components names must be given
            error('IntrepidTwilight:executive:build:NotAllObjectChosen',...
                'Component ''%s'' could not be built because not all object choices are set',...
                componentName);
        end
        
        
    end
    
end