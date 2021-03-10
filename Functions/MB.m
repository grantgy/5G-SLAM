classdef MB
    
    properties
        p_exist;
        Gaussian_mixture;
        log_weight;
        birth_time
    end
    
    methods
        function obj = MB(a,b,c,d)
            obj.p_exist = a;
            obj.Gaussian_mixture = b;
            obj.log_weight = c;
            obj.birth_time = d;
        end
    end
end
