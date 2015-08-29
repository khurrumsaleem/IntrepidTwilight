function solver = Solver(~)
    
    solver = IntrepidTwilight.executive.Component();
    solver = solver.changeID(solver,'solver','solver');
    
    %   Base dependencies
    solver.dependencies = {'residual','preconditioner'};
    
    %   Identity solve
    solver.solve = @(v) v;
    
end