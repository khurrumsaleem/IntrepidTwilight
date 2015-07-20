function closure = BaseClosure(type)

    closure      = IntrepidTwilight.executive.Registry();
    closure.type = type;
    closure.is   = @(s) strcmpi(s,type);

end