function pc = Preconditioner(residual,type)
    
    r    = residual ;
    drdq = []       ;
    
    
    switch(lower(type))

        case('full-stagnant')

            pc = struct(...
                'apply'      , @(q) applyFullJacobian(q)        ,...
                'initialize' , @(q) initializeFullJacobian(q)   ,...
                'update'     , @(q) []                          ,...
                'get'        , @() get()                        );

        case('block-jacobi')

            pc = struct(...
                'apply'      , @(q) applyBlockJacobi(q)     ,...
                'initialize' , @(q) updateBlockJacobi(q)    ,...
                'update'     , @(q) updateBlockJacobi(q)    ,...
                'get'        , @() get()                    );

        case('none')

            pc = struct(...
                'apply'      , @(q) q    ,...
                'initialize' , @(q) q    ,...
                'update'     , @(q) []   ,...
                'get'        , @()  []   );

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
    function j = get()
        j = drdq;
    end



    % ======================================================================= %
    %                         Staganat Full functions                         %
    % ======================================================================= %
    function u = applyFullJacobian(q)
        u = drdq \ q;
    end
    function [] = initializeFullJacobian(q)
        drdq = r.jacobian(q);
    end
    
    
end