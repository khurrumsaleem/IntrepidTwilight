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
    component = 'state'                                                     ;
    config    = components.(component).get()                                ;
    name      = config.name                                                 ;
	state     = IntrepidTwilight.executive.build('executive',name,config)   ;
    state     = state.changeID(state,'HEMFull')                             ;
    
    
    %   Array used to store  previously calculated quantities
    %   for use in linearity estimation
    tstore    = [];
    qstore    = [];
    Dqstore   = [];
    isDynamic = [];
    
    
    %   Prepare
    state.prepare = @(varargin) prepare(varargin{:});
    isNotPrepared = true;
    warehouse     = [];
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
            isDynamic     = sd.get('isDynamic')     ;
            isNotPrepared = false                   ;
            
            %   Pull all state-set information
            warehouse = state.get();

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
            
            %   Offer a hook before anything changes
            warehouse.hook.preupdate(q,t,dt);
            
            
            %   Update residual
            r.update(q,t,dt);

            
            % -----------------------------------------------------------------
            %   Perform a Linearly Implicit Euler step in an attempt to 
            %   find a good nonlinear search direction.
            q  = IntrepidTwilight.ConvenientMeans.linearlyImplicitEuler(...
                @(q,t) sd.rhs(q,t),q,linspace(0,dt,2));
            %
            %   Adjust for dynamism
            dq = (q(:,1)-q(:,end)).*isDynamic;
            %
            %   Scale the LIE step since it is not nonlinearly stable.
            if not(any(isnan(dq))) && not(any(abs(imag(dq)) > eps()))
                dq = r.guard.step(q(:,1),dq);
                while any(abs(dq)>1E-12) && any(isnan(r.value(q(:,1) - dq)))
                    dq = 0.5 * dq;
                end
                q = q(:,1) - dq;
            else
                q = q(:,1);
            end
            % -----------------------------------------------------------------
            
            %   Solve
            [q,stats] = solver.solve(q);
            
            %   Update data while accouting for a fallback
            tstore  = [t,tstore(1)]                     ;
            qstore  = [q,qstore(:,1)]                   ;
            Dqstore = [sd.rhs(q,t+dt),Dqstore(:,1)]     ;
            Dq      = sd.makeDimensional(Dqstore(:,1))  ;
            q       = sd.makeDimensional(qstore(:,1))   ;


            %   Offer a hook after the potentially new q is found
            warehouse.hook.postupdate(q,t,dt,stats);

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
