function [f,varargout] = SemidiscreteUpwind(q,s)
    
    [frho ,s] = SemidiscreteUpwindMass(q,s);
    [frhoe,s] = SemidiscreteUpwindEnergy(q,s);
    [frhov,s] = SemidiscreteUpwindMomentum(q,s);
    
    
    f = [frho;frhoe;frhov];
    
    if (nargout > 1)
        varargout{1} = s;
    end

end