classdef global_hypothesis
    
    properties
        weight;
        look_up_table;
    end
    
    methods
        function obj = global_hypothesis(a,b)
            obj.weight = a;
            obj.look_up_table = b;
        end
    end
end
