function hem = HEM()
    
    %   Inherit
    hem = IntrepidTwilight.executive.Simulation();
    hem = hem.changeID(hem,'simulation','HEM');
    
    %   Set defaults
    hem.choose('model'              , 'AdamantWave'        , 'BasicModel'       );
    hem.choose('semidiscretization' , 'AdamantWave'        , 'Quasi2DUpwind'    );
    hem.choose('timediscretization' , 'TradientStride'     , 'ImplicitEuler'    );
    hem.choose('residual'           , 'executive'          , 'Residual'         );
    hem.choose('preconditioner'     , 'executive'          , 'Preconditioner'   );
    hem.choose('solver'             , 'TenaciousReduction' , 'JFNK'             );
    hem.choose('state'              , 'executive'          , 'State'            );
    hem.choose('evolver'            , 'executive'          , 'Evolver'          );
    
    evolver = [];
    hem.run = @() run();
    function [] = run()
        
        %   Build everything (done now to ensure all configuration is done)
        %         evolver = buildNestedEvolver(hem);
        %         evolver = buildSerialEvolver(hem);
        evolver = buildFullEvolver(hem);
        
        
        %   Run it.
        evolver.run();
        
    end
    
    hem.results = @() results();
    function varargout = results()
        varargout = cell(1,nargout);
        [varargout{:}] = evolver.getData();
    end
    
    
end


function evolver = buildFullEvolver(hem)
    
    %   Get all simulation-set data
    components = hem.get();
    
    
    
    % ======================================================================= %
    %                          Non-Solver Components                          %
    % ======================================================================= %
    
    % ----------------------------------------- %
    %                     Model                 %
    % ----------------------------------------- %
    component = 'model';
    config    = components.(component).get();
    name      = config.name;
    model     = IntrepidTwilight.executive.build('AdamantWave',name,config);
    
    
    
    % ----------------------------------------- %
    %            Semidiscretizations            %
    % ----------------------------------------- %
    
    %   Canonical semidiscretization
    component = 'semidiscretization';
    config    = components.(component).get();
    name      = config.name;
    sd = IntrepidTwilight.executive.build('AdamantWave',name,config);
    sd.bind(model);
    
    
    
    % ----------------------------------------- %
    %            Timediscretizations            %
    % ----------------------------------------- %
    component = 'timediscretization'        ;
    config    = components.(component).get();
    name      = config.name                 ;
    ts        = IntrepidTwilight.executive.build('TransientStride',...
        name,config);
    ts.bind(sd);
    
    
    
    % ----------------------------------------- %
    %                  Residuals                %
    % ----------------------------------------- %
    component = 'residual'                  ;
    config    = components.(component).get();
    name      = config.name                 ;
    r = IntrepidTwilight.executive.build('executive',...
        name,config);
    r.bind(ts);
    
    
    % ----------------------------------------- %
    %               Preconditioners             %
    % ----------------------------------------- %
    component = 'preconditioner'            ;
    config    = components.(component).get();
    name      = config.name                 ;
    pc = IntrepidTwilight.executive.build('executive',...
        name,config);
    pc.bind(r);
    
    
    
    
    
    % ======================================================================= %
    %                                 Solvers                                 %
    % ======================================================================= %
    
    %   Build
    component = 'solver'                    ;
    config    = components.(component).get();
    name      = config.name                 ;
    solver = IntrepidTwilight.executive.build('TenaciousReduction',...
        name,config);
    solver.bind(r,pc);

    
    
    % ======================================================================= %
    %                                 State                                   %
    % ======================================================================= %
    state = IntrepidTwilight.executive.State()  ;
    state = state.changeID(state,'HEMFull')     ;
    
    
    %   Array used to store  previously calculated quantities
    %   for use in linearity estimation
    tstore  = [];
    qstore  = [];
    Dqstore = [];
    
    
    %   Prepare
    state.prepare = @(varargin) prepare(varargin{:});
    isNotPrepared = true;
    function [] = prepare(q,t,varargin)
        if isNotPrepared
            
            %   Prepare cascade
            r.prepare(q,t,varargin{:});
            pc.prepare();
            solver.prepare(q,t,varargin{:});
            
            %   Store IC and block preparation
            tstore        = t                       ;
            qstore        = sd.makeDimensionless(q) ;
            Dqstore       = sd.rhs(qstore,t)        ;
            isNotPrepared = false                   ;
            
        end
    end
    
    
    % ----------------------------------------- %
    %                   Update                  %
    % ----------------------------------------- %
    state.update = @(q,t,dt) update(q,t,dt);
    function [q,Dq,stats] = update(q,t,dt)
        
        %   Split data
        q = sd.makeDimensionless(q);
        
        if (abs(dt) > 0)
            %   Update residual
            r.update(q,t,dt);
            %
            %   Perform a simple Euler step to knock solution off of current state
            if (dt < 1E-6)
%                 q = q + dt*sd.rhs(q,t);
            else
%                 q = q + 1E-4*sd.rhs(q,t);
            end
            

            %   Solve
            [q,stats] = solver.solve(q);

            
            %   Update data while accouting for a fallback
            tstore  = [t,tstore(1)]                     ;
            qstore  = [q,qstore(:,1)]                   ;
            Dqstore = [sd.rhs(q,t+dt),Dqstore(:,1)]     ;
            Dq      = sd.makeDimensional(Dqstore(:,1))  ;
            q       = sd.makeDimensional(qstore(:,1))   ;

        else
            Dq      = sd.makeDimensional(sd.rhs(q,t))   ;
            q       = sd.makeDimensional(q)             ;
        end
        

        
        
        %   Use linearity estimate as error estimate for the Implicit
        %   Euler solution
        stats.absoluteErrorEstimate = dt*norm(diff(Dqstore,[],2))               ;
        stats.relativeErrorEstimate = stats.absoluteErrorEstimate/norm(qstore)  ;
        
    end
    
    
    % ======================================================= %
    %                         Evolver                         %
    % ======================================================= %
    component = 'evolver'                   ;
    config    = components.(component).get();
    name      = config.name                 ;
    evolver   = IntrepidTwilight.executive.build('executive',...
        name,config);
    evolver.bind(state);

    
    
end
