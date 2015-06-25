function semid = makeSemidiscretization(model)

    semidName = model.semidiscretization.name;
    modules   = IntrepidTwilight.executive.moduleRegistry('physics');
    
    for k = 1:length(modules)
        module = modules{k};
        if which(['IntrepidTwilight.',module,'.Semidiscretizations.registry'])
            semids = IntrepidTwilight.(module).Semidiscretizations.registry();
            
            if isfield(semids,semidName)
                break;
            end
        end
    end
    
    semid = IntrepidTwilight.(module).Semidiscretizations.(semidName)(model);

end