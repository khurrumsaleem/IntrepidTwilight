function pc = Preconditioner(problem)
    
    
    switch(lower(problem.solver.preconditioner.type))

        case('full-stagnant')
            
            %   Define preconditioner closure
            pc = @(dt) struct(...
                'apply'      , @(q) applyFullJacobian(q)        ,...
                'initialize' , @(q) initializeFullJacobian(q,dt),...
                'update'     , @(q) []                          ,...
                'get'        , @() get());
        

        case('block-jacobi')
            
            
            %   Build economy identity matrix
            blockSize = problem.solver.preconditioner.blockSize;
            rEye = arrayfun(@(e) eye(e),blockSize,'UniformOutput',false);

            %   Get initial dfdq
            dfdq = ...
                problem.semidiscretization.closure.blockDiagonalJacobian(problem.initialState.q0);
            
            %   Initialize drdq
            drdq = rEye;
            
            %   Define preconditioner closure
            pc = @(dt) struct(...
                'apply'      , @(q) applyBlockJacobi(q),...
                'initialize' , @(q) updateBlockJacobi(q,dt),...
                'update'     , @(q) updateBlockJacobi(q,dt),...
                'get'        , @() get());
            
        case('none')
            pc = @(dum) struct('apply',@(x)x,'update',@(dum)[],'get',@()[]);
            
    end
    
    
    
    % ======================================================================= %
    %                         Block-Jacobi functions                          %
    % ======================================================================= %
    function u = applyBlockJacobi(q)
        u = IntrepidTwilight.TenaciousReduction.blockDiagonalEconomy(drdq,q,blockSize);
    end
    function [] = updateBlockJacobi(q,dt)
        dfdq = problem.semidiscretization.closure.blockDiagonalJacobian(q)      ;
        drdq = cellfun(@(I,dF) I - dt*dF , rEye , dfdq,'UniformOutput',false)   ;
    end
    function j = get()
        j = drdq;
    end



    % ======================================================================= %
    %                         Block-Jacobi functions                          %
    % ======================================================================= %
    function u = applyFullJacobian(q)
        u = drdq \ q;
    end
    function [] = initializeFullJacobian(q,dt)
        dfdq = IntrepidTwilight.ConvenientMeans.numericalJacobianFull(...
            problem.semidiscretization.closure.rhs,q);
        drdq = eye(size(dfdq)) - dt*dfdq                                ;
    end
    
    
end