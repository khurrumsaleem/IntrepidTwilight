function state = State(config)
    
    %   Inherit
    state = IntrepidTwilight.executive.Component();
    state = state.changeID(state,'state','state');
    
    %    Parameters
    if (nargin >= 1) && not(isempty(config))
        state.set(config);
    end



    %   Binder
    residual = [];
    solver   = [];
    state.set('dependencies',{'residual','solver'});
    state.bind = @(varargin) bind(varargin);
    function [] = bind(objects)
        for k = 1:length(objects)
            bind_(objects{k});
        end
    end
    function [] = bind_(object)
        switch(lower(object.type))
            case('residual')
                residual = object;
            case('solver')
                solver = object;
        end
    end



    %   Methods
    state.update = @(q,t,dt) update(q,t,dt);
    function [q,t,stats] = update(q,t,dt)
        
        %   Update residual
        residual.update(q,t,dt);
            
        %   Solve
        [q,stats] = solver.solve(q);
    end
    
    
    %   Prepare
    state.prepare = @(varargin) prepare(varargin{:});
    function [] = prepare(q,t,varargin)
        %   Preparation cascade
        solver.prepare(q,t,varargin{:});
    end
    
    
end