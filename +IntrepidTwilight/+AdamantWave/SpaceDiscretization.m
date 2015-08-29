function sd = SpaceDiscretization(~)

    %   Inherit and overwrite ID information
    sd = IntrepidTwilight.executive.Component() ;
    sd = sd.changeID(sd,'spaceDiscretization','spaceDiscretization');
    
    
%     %   Binder
%     model   = [] ; 
%     sd.bind = @(m) bind(m)  ;
%     function [] = bind(m)
%         model = m;
%     end


end