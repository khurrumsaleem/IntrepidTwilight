function [components,modules,names] = defaultCMNLists()

    components  = { 'model','spacediscretization','timediscretization',...
                    'residual','preconditioner','solver','evolver'};
    modules     = { '','','TransientStride','executive','executive',...
                    'TenaciousReduction','executive'};
    names       = { '','','ImplicitEuler','Residual','Preconditioner',...
                    'JFNK','Evolver'};

end