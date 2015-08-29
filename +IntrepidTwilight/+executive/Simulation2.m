function simulation = Simulation2()
    
    %   Set-Up
    simulation = IntrepidTwilight.executive.Object();
    simulation = simulation.changeID(simulation,'simulation','simulation');
    simulation.set('bindAtBuild',true);
    
    
    %   Default values
    [components,modules,names] = IntrepidTwilight.executive.defaultCMNLists();
    
    
    for k = 1:numel(components)
        component = components{k};

        %   Make setter/getter
        simulation.(component).set = @(key,value) simulation.set([component,'.',key],value) ;
        simulation.(component).get = @(varargin)  simulation.get([component,varargin{:}])   ;
        
        %   Set default values
        simulation.(component).set('module',modules{k});
        simulation.(component).set('name'  ,names{k});
    end
    
    
    %   Define build
    simulation.build = @() build();
    function [] = build()
        componentFun(@(s,c) buildLocal(s,c),simulation,components);
    end


    
    
    %   Define simple evolver
    simulation.evolve = @() evolve();
    function [] = evolve()
        componentFun(@(s,c)prepare(s,c) , simulation , components);
    end
    
end

function [] = componentFun(fun,simulation,componentNames)
    for k = 1:numel(componentNames)
        fun(simulation,componentNames{k});
    end
end


function [] = buildLocal(simulation,componentName)
    
    component = simulation.get(componentName)   ;
    module    = component.module                ;
    name      = component.name                  ;


    if not(isempty(module) || isempty(name))

        %   (Attempt to) Build
        object = IntrepidTwilight.executive.build(module,name);


        %   (Attempt to) Inject dependencies
        if  store.bindAtBuild               && ...
            isfield(object,'dependencies')  && ...
            not(isempty(object.dependencies))

            dependencies = object.dependencies;
            for m = 1:numel(dependencies)
                object.bind(simulation.get(dependencies{m}));
            end

        end

        %   Set in store
        simulation.set(component,object);

    else

        error('IntrepidTwilight:executive:Simulation:objectNotFound',...
            ['The component ''%s'' could not be built because the object''s module or ',...
            'name was unspecified.']);

    end

end

function [] = prepare(simulation,componentName)
    component = simulation.(componentName);
    component.prepare();
end


