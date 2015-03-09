function c = blah()
   
    a = 1;
    
    c.add = @add;
    c.get = @get;
    
    function [] = add(x)
        a = a + x;
    end
    function out = get()
        out = a;
    end
    
    
    
end