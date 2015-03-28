function pc = Preconditioner(problem)
    
    
    switch(lower(problem.solver.preconditioner.type))

        case('block-jacobi')
            
            
            blockSize = problem.solver.preconditioner.blockSize;
            nColumns  = max(blockSize)                  ;
            nBlocks   = numel(blockSize)                ;

            
            %   Build economy identity matrix
            s         = blockSize(1)                                ;
            J         = 1:s                                         ;
            I         = J                                           ;
            rEye      = zeros(problem.miscellaneous.nEq,nColumns)   ;
            rEye(I,J) = eye(blockSize(1))                           ;
            for k = 2:nBlocks
                J         = 1:blockSize(k)      ;
                I         = s + J               ;
                rEye(I,J) = eye(blockSize(k))   ;
                s         = I(end)              ;
            end


            %   Get initial dfdq
            dfdq = ...
                problem.semidiscretization.closure.blockDiagonalJacobian(problem.initialState.q0);
            
            %   Initialize drdq
            drdq = rEye;
            
            %   Define preconditioner closure
            pc = @(dt) struct('apply',@(q) applyBlockJacobi(q),'update',@(q) updateBlockJacobi(q,dt));
            
        case('none')
            pc.apply  = @(q) q  ; 
            pc.update = @(q) [] ;
            
    end
    
    
    
    % ======================================================================= %
    %                         Block-Jacobi functions                          %
    % ======================================================================= %
    function u = applyBlockJacobi(q)
        u = IntrepidTwilight.TenaciousReduction.blockDiagonalEconomy(drdq,q,blockSize);
    end
    function [] = updateBlockJacobi(q,dt)
        dfdq = problem.semidiscretization.closure.blockDiagonalJacobian(q);
        drdq = rEye - dt*dfdq;
    end
    
    
end