function [Newmap] = mapreconstruction(local_hyp,weight,Gat_treshold,Prun_treshold,map_BS,global_hyp)

    VA_pos = []; SP_pos = [];
    VA_P = []; SP_P = [];
    VA_weight = []; SP_weight = [];
    
    for i = 1:size(local_hyp,2)
        VA_weight_eachpart =[];SP_weight_eachpart =[];
        
        for j = 2:size(local_hyp{1,i},2)
        
            for k = 1:size(local_hyp{1,i}{1,j}.p_exist,2)
                w = sum(global_hyp{1,i}.weight(1,find(global_hyp{1,i}.look_up_table(:,j))==k));
                VA_pos = [VA_pos,local_hyp{1,i}{1,j}.Gaussian_mixture(1,k).mean(:,1)];
                VA_P = cat(3,VA_P,local_hyp{1,i}{1,j}.Gaussian_mixture(1,k).covariance(:,:,1));
                VA_weight_eachpart = [VA_weight_eachpart,local_hyp{1,i}{1,j}.p_exist(1,k)*local_hyp{1,i}{1,j}.Gaussian_mixture(1,k).weight(1,1)*w];
                
                if length(local_hyp{1,i}{1,j}.Gaussian_mixture(1,k).weight) > 1
                    SP_pos = [SP_pos,local_hyp{1,i}{1,j}.Gaussian_mixture(1,k).mean(:,2)];
                    SP_P = cat(3,SP_P,local_hyp{1,i}{1,j}.Gaussian_mixture(1,k).covariance(:,:,2));
                    SP_weight_eachpart = [SP_weight_eachpart,local_hyp{1,i}{1,j}.p_exist(1,k)*local_hyp{1,i}{1,j}.Gaussian_mixture(1,k).weight(1,2)*w];
                end
            end
        end
        
        VA_weight = [VA_weight,VA_weight_eachpart*weight(1,i)];
        SP_weight = [SP_weight,SP_weight_eachpart*weight(1,i)];
    end
        
    [Newmap_VA] = map_pruning(VA_pos, VA_P, VA_weight,Gat_treshold,Prun_treshold);
    [Newmap_SP] = map_pruning(SP_pos, SP_P, SP_weight,Gat_treshold,Prun_treshold);
    Newmap = MAP(map_BS,Newmap_VA,Newmap_SP); 
end

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
                index_temp =[index_temp,i];
            end
        end
        
        weight_temp = weight(1,index_temp);
        position_temp = position(:,index_temp);
        covariance_temp = covariance(:,:,index_temp);        
        weight_temp_sum = sum(weight_temp);
        New_weight = [New_weight,weight_temp_sum];      
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
