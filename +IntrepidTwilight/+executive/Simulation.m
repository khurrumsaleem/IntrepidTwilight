function simulation = Simulation(moduleName)
    
    %   Inherit
    simulation = IntrepidTwilight.executive.Object();
    simulation = simulation.changeID(simulation,'simulation','simulation');
    
    if (nargin >= 1) && not(isempty(moduleName))
        setBuildModuleName(moduleName);
    end
    
    
    %   Build order
    simulation.set('buildOrder',...
        {'model','spacediscretization','timediscretization','residual',...
        'preconditioner','solver','state'});
    
    
    %   Storage struct for component configurations that can be specified
    %   after a component is chosen.
    components     = struct();
    
    
    %   Public methods
    simulation.build  = @(varargin) build(varargin{:})  ;
    simulation.evolve = @() evolve()                    ;
    simulation.run    = @() evolve()                    ;
    simulation.data   = @() getData()                   ;
    
    %   Late module name bind
    simulation.setBuildModuleName = @(component,name) setBuildModuleName(component,name);
    
    
    
    
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
    simulation.bind = @(component) bind(component);
    function [] = bind(component)
        if isstruct(component)
            components.(component.type)        = component  ;
            components.(component.type).config = ...
                IntrepidTwilight.executive.Object()         ;
        end
    end


    %   Choose scheme
    simulation.choose = @(component,module,name) choose(component,module,name)  ;
    function [] = choose(component,module,name)
        components.(component) = IntrepidTwilight.executive.Object();
        components.(component).set('module',module)                 ;
        components.(component).set('name'  ,name)                   ;
    end


    %   Overloaded getter
    simulation.get = @(varargin) get(varargin{:});
    function value = get(component,varargin)
        if (nargin > 1) && isfield(components,component)
            value = components.(component).get(varargin{:});
        elseif (nargin == 1)
            value = components.(component);
        elseif (nargin == 0)
            value = components;
        end
    end

    
    %   Overloaded setter
    simulation.set = @(component,varargin) set(component,varargin{:});
    function [] = set(component,varargin)
        if isfield(components,component)
            components.(component).set(varargin{:});
        end
    end





    % =================================================== %
    %                    Build System                     %
    % =================================================== %
    
    function [] = build(bindAtBuild)
        
        if (nargin < 1) || isempty(bindAtBuild)
            bindAtBuild = true;
        end
        
        buildOrder = simulation.get('buildOrder');
        for component = buildOrder

            if isfield(configurations,component)
                
                %   Build component
                module = components.(component).get([component,'.module']);
                name   = components.(component).get([component,'.name']);
                config = components.(component).get();
                components.(component) = ...
                    IntrepidTwilight.executive.build(module,name,config);


                %   Bind if requested
                if bindAtBuild
                    dependencies = components.(component).get('dependencies');
                    if not(isempty(dependencies))
                        for dependency = 1:numel(dependencies)
                            if isfield(components,dependency)
                                components.(component).bind(...
                                    components.(dependency));
                            end
                        end
                    end
                end

            end
        end
    end
    
end