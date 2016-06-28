function [xs,stats] = serialSolver(xs,solvers,continueSweep)
    
    %   Default postFlags predicate
    if (nargin < 3) || isempty(continueSweep) || not(isa(continueSweep,'function_handle'))
        continueSweep = @(xs,stats,postFlags)           ...
            not(cellfun(@(c) c.converged,stats))    &&  ...
            any(strcmpi('NotDone',postFlags))       ;
    end


    %   Loop set-up
    n              = numel(xs)  ;
    stats{1,n}     = []         ;
    postFlags{1,n} = []         ;
    notDone        = true       ;
    iterationMax   = 100        ;
    iteration      = 0          ;


    %   Iterative
    while notDone
        
        %   Serial Sweep
        for k = 1:n
            [xs{k},stats{k},postFlags{k}] = solvers{k}.solve(xs{k});
        end
        
        %   Post sweep stuff
        iteration    = iteration + 1                                    ;
        belowIterMax = iteration < iterationMax                         ;
        notDone      = belowIterMax & continueSweep(xs,stats,postFlags) ;
        
    end


end