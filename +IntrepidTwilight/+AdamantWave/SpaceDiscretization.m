function sd = SpaceDiscretization(~)

    %   Inherit and overwrite ID information
    sd = IntrepidTwilight.executive.Component() ;
    sd = sd.changeID(sd,'spaceDiscretization','spaceDiscretization');
    sd.set('dependencies',{'model'});

    sd.rhs     = @(q,t) [];
    sd.getAll  = @()    [];
    sd.prepare = @(varargin)  [];

end