function pc = Preconditioner(residual,kind)
    
    %   Closure variables
    r    = []   ;
    drdq = []   ;
    
    
    %   Instantiate struct
    pc.is   = @(s) strcmpi(s,'preconditioner')  ;
    pc.get  = @()  get()                        ;
    pc.set  = @(type,object) set(type,object)   ;
    
    
    %   Bind at construction
    if (nargin >= 1)
        r = residual;
    end
    if (nargin >= 2) && not(isempty(kind))
        bindMethods();
    end
    
    
    
    function [] = bindMethods()
        
        if r.is('residual')
            switch(lower(kind))
                
                case('full-stagnant')
                    pc.apply      = @(q) applyFullJacobian(q)        ;
                    pc.initialize = @(q) initializeFullJacobian(q)   ;
                    pc.update     = @(q) []                          ;
                    
                    
                case('block-jacobi')
                    pc.apply      = @(q) applyBlockJacobi(q)     ;
                    pc.initialize = @(q) updateBlockJacobi(q)    ;
                    pc.update     = @(q) updateBlockJacobi(q)    ;
                    
                    
                case('none')
                    pc.apply      = @(q) q   ;
                    pc.initialize = @(q) []  ;
                    pc.update     = @(q) []  ;
                    
            end
        end
    end
    
    % ======================================================================= %
    %                                Re-binders                               %
    % ======================================================================= %
    function [] = set(type,object)
        switch(lower(type))
            case('residual')
                if object.is('residual')
                    r = object;
                end

            case('kind')
                kind = object;
                bindMethods();
        end
    end
    function j = get()
        j = drdq;
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