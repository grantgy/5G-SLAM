classdef GM
    
    properties
        weight;
        mean;
        covariance;
    end
    
    methods
        function obj = GM(a,b,c)
            obj.weight = a;
            obj.mean = b;
            obj.covariance = c;
        end
    end
end