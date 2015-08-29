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
    pc.dependencies = {'residual'}          ;
    pc.bind         = @(object) bind(object);
    %
    function [] = bind(object)
        if isstruct(object) && object.is('residual')
            r = object;
        end
    end



    %   Prepare for solution
    pc.prepare = @() prepare();
    %
    function [] = prepare()
        if r.is('residual')
            bindMethods();
        else
            error('IntrepidTwilight:executive:Preconditioner:noResidual',...
                'Preconditioner does not have a residual bound to it.');
        end
    end


    function [] = bindMethods()

        switch(lower(pc.get('kind')))
                
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