function component = Component()
    
    component = IntrepidTwilight.executive.Object();
    component = component.changeID(component,'component','component');
    
    %   Dependency store
    component.dependencies = {};
    
    %   Prepare handle
    component.prepare = @() [];
    
end