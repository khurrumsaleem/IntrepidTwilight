function pc = Preconditioner(problem)
    
    
    switch(lower(problem.solver.preconditioner.type))

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
                'apply',@(q) applyBlockJacobi(q),...
                'update',@(q) updateBlockJacobi(q,dt),...
                'get',@() get());
            
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
    
    
end