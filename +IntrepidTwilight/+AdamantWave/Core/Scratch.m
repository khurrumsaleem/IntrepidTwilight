function xSol = Scratch()
    
    Update = @(x,Mask) SquareRoot(x,Mask);
    xSol   = NewtonUpdater(Update,[2;5;12;40056]+rand(4,1),1E-17,100);

end

function [dx,R] = SquareRoot(x,~)
    
    
    f    = x.^(1/3);
    dfdx = (1/3).*x.^(-2/3);
    dx   = f./dfdx;
    R    = abs(real(dx)) + imag(f)*1i;
    
    
end