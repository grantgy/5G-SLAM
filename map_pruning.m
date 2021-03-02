function [New_map] = map_pruning(position,covariance,weight,Gat_treshold,Prun_treshold)

    New_position=[];
    New_covariance=[];
    New_weight=[];
    
    while ~isempty(weight)
        index = find(weight == max(weight));
        index = index(1,1);
        
        index_temp = [];
        for i = 1:length(weight)
            if norm(position(:,index)-position(:,i)) <  5
%             if (position(:,index)-position(:,i))'/covariance(:,:,i)*(position(:,index)-position(:,i)) <  Gat_treshold
                index_temp =[index_temp,i];
            end
        end
        weight_temp = weight(1,index_temp);
        position_temp = position(:,index_temp);
        covariance_temp = covariance(:,:,index_temp);
        
        weight_temp_sum = sum(weight_temp);
        New_weight = [New_weight,weight_temp_sum];
        %New_weight = [New_weight,weight(index)];
        
        weight_temp = weight_temp./weight_temp_sum;
        
        [x_new, P_new] = merge(GM(weight_temp,position_temp,covariance_temp));

        New_position=[New_position,x_new];
        New_covariance=cat(3,New_covariance,P_new);
        
        index_new = setdiff((1:length(weight)),index_temp);
        
        weight = weight(1,index_new);
        position = position(:,index_new);
        covariance = covariance(:,:,index_new);
        
    end
    
    if ~isempty(New_weight)
        
        index1 = find(New_weight > Prun_treshold);
        New_position=New_position(:,index1);
        New_covariance=New_covariance(:,:,index1);
        New_weight=New_weight(1,index1);

    end
    New_map = GM(New_weight,New_position,New_covariance);
end
