function pc = Preconditioner(config)
    
    %   Instantiate struct
    pc = IntrepidTwilight.executive.Component();
    pc = pc.changeID(pc,'preconditioner','preconditioner');
    
    
    
    %   Set defaults
    pc.set('kind','none');
    %
    %   Overwrite defaults at construction
    if (nargin >= 1)
        pc.set(config);
    end
    
    
    
    %   Dependencies and Binder
    r               = []                    ;
    drdq            = []                    ;
    pc.set('dependencies',{'residual'})     ;
    pc.bind         = @(object) bind(object);
    %
    function [] = bind(object)
        if isstruct(object) && object.is('residual')
            r = object;
        end
    end
    
    
    
    %   Prepare for solution
    pc.prepare    = @(varargin) prepare(varargin{:});
    isNotPrepared = true;
    function [] = prepare(varargin)
        if isNotPrepared
            if r.is('residual')
                bindMethods();
            else
                error('IntrepidTwilight:executive:Preconditioner:noResidual',...
                    'Preconditioner does not have a residual bound to it.');
            end
            isNotPrepared = false;
        end
    end
    
    
    
    %   Define dynamic method binding interface
    apply_        = [];
    initialize_   = [];
    update_       = [];
    pc.apply      = @(q) apply(q);
    pc.initialize = @(q) initialize(q);
    pc.update     = @(q) update(q);
    function out = apply(q)
        out = apply_(q);
    end
    function [] = initialize(q)
        initialize_(q);
    end
    function [] = update(q)
        update_(q);
    end
    
    
    
    %   Define run-time binder
    function [] = bindMethods()
        
        switch(lower(pc.get('kind')))
            case('full-stagnant')
                apply_      = @(q) applyFullJacobian(q)        ;
                initialize_ = @(q) initializeFullJacobian(q)   ;
                update_     = @(q) []                          ;
                
            case('block-jacobi')
                apply_      = @(q) applyBlockJacobi(q)     ;
                initialize_ = @(q) updateBlockJacobi(q)    ;
                update_     = @(q) updateBlockJacobi(q)    ;
                
            case('none')
                apply_      = @(q) q   ;
                initialize_ = @(q) []  ;
                update_     = @(q) []  ;
        end
        
    end
    
    
    % ======================================================================= %
    %                         Block-Jacobi functions                          %
    % ======================================================================= %
    function u = applyBlockJacobi(q)
        u = IntrepidTwilight.TenaciousReduction.blockDiagonalEconomy(drdq,q);
    end
    function [] = updateBlockJacobi(q)
        drdq = r.blockDiagonalJacobian(q);
    end
    
    
    
    % ======================================================================= %
    %                         Stagnanat Full functions                        %
    % ======================================================================= %
    function u = applyFullJacobian(q)
        u = drdq \ q;
    end
    function [] = initializeFullJacobian(q)
        drdq = r.jacobian(q);
    end
    
    
end