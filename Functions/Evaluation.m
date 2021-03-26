function [error_location,error_heading,error_bias,GOSPA_VA,GOSPA_SP] =Evaluation(state,real_state,map,BS,VA,SP,p, c, alpha,T_VA,T_SP)
    error_location =  sqrt(sum((state(1:2,:) - real_state(1:2,:)).^2,1));
    error_heading = abs(state(3,:) - real_state(3,:));
    
    for i = 1:size(error_heading,2)
        while error_heading(1,i) > 2*pi
            error_heading(1,i) =error_heading(1,i)-2*pi; 
        end
        if error_heading(1,i) > pi
            error_heading(1,i) = 2*pi - error_heading(1,i); 
        end
    end
    
    error_bias = abs(state(4,:) - real_state(4,:));
    GOSPA_VA = zeros(1,size(map,2));
    GOSPA_SP = zeros(1,size(map,2));
    
    for i = 2 : size(map,2)
        if ~isempty(map{1,i}.VA.mean) && ~isempty(VA) 
            map{1,i}.VA.mean =  map{1,i}.VA.mean(:,map{1,i}.VA.weight >= T_VA);
            [GOSPA_VA(1,i),~, ~] = GOSPA(map{1,i}.VA.mean, VA, p, c, alpha);
        elseif isempty(map{1,i}.VA.mean) && ~isempty(VA) 
            GOSPA_VA(1,i)= GOSPA(double.empty(3,0), VA, p, c, alpha);
        elseif ~isempty(map{1,i}.VA.mean) && isempty(VA)   
            map{1,i}.VA.mean =  map{1,i}.VA.mean(:,map{1,i}.VA.weight >= T_VA);
            GOSPA_VA(1,i)= GOSPA(map{1,i}.VA.mean,double.empty(3,0), p, c, alpha);
        else
            GOSPA_VA(1,i)= GOSPA(double.empty(3,0), double.empty(3,0), p, c, alpha);
        end
    end
    
    for i = 2:size(map,2)
        if ~isempty(map{1,i}.SP.mean) && ~isempty(SP)
            map{1,i}.SP.mean =  map{1,i}.SP.mean(:,map{1,i}.SP.weight >= T_SP);
            [GOSPA_SP(1,i),~, ~] = GOSPA(map{1,i}.SP.mean, SP, p, c, alpha);
        elseif isempty(map{1,i}.SP.mean) &&  ~isempty(SP)
            GOSPA_SP(1,i) = GOSPA(double.empty(3,0), SP, p, c, alpha);
        elseif ~isempty(map{1,i}.SP.mean) &&  isempty(SP)
            map{1,i}.SP.mean =  map{1,i}.SP.mean(:,map{1,i}.SP.weight >= T_SP);
            GOSPA_SP(1,i) = GOSPA(map{1,i}.SP.mean,double.empty(3,0), p, c, alpha);
        else
            GOSPA_SP(1,i)= GOSPA(double.empty(3,0), double.empty(3,0), p, c, alpha);
        end
    end
end
