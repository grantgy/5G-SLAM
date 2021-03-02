classdef MAP
    
    properties
        BS;
        VA;
        SP;
    end
    
    methods
        function obj = MAP(a,b,c)
            obj.BS = a;
            obj.VA = b;
            obj.SP = c;
        end
    end
end