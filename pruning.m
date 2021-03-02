function [PPP_pruned,local_hyp_pruned,global_hyp_pruned] = pruning(PPP,local_hyp,global_hyp,T_pruning,T_pruningPois,Nhyp_max,existence_threshold)
    PPP_pruned = PPP;
    local_hyp_pruned = local_hyp;
    global_hyp_pruned = global_hyp;
    
    % prune PPP
    for i = 1:length(PPP)
        index_remove = PPP(i,1).weight < T_pruningPois;
        PPP_pruned(i,1).weight(index_remove) = [];
        PPP_pruned(i,1).pro(index_remove) = [];
%         PPP_pruned(i,1).mean(:,index_remove) = [];
%         PPP_pruned(i,1).covariance(:,:,index_remove) = [];
        
    end
    % prune global_hyp_pruned
    index_remove = global_hyp.weight < T_pruning;
    global_hyp_pruned.weight(index_remove) = [];
    global_hyp_pruned.look_up_table(index_remove,:) = [];
    
    if length(global_hyp_pruned.weight) > Nhyp_max
        [~,index_sorted] = sort(global_hyp_pruned.weight,'descend');
        global_hyp_pruned.weight= global_hyp_pruned.weight(index_sorted(1:Nhyp_max));
        global_hyp_pruned.look_up_table = global_hyp_pruned.look_up_table(index_sorted(1:Nhyp_max),:);
    end
    global_hyp_pruned.weight = global_hyp_pruned.weight./sum(global_hyp_pruned.weight);
    
    % remove object
    index_remove = find (sum(global_hyp_pruned.look_up_table,1)==0);
    local_hyp_pruned(:,index_remove) = [];
    global_hyp_pruned.look_up_table(:,index_remove) = [];
    
    % remove local hypothesis according to global hypothesis
    for i = 1:length(local_hyp_pruned)
        index_valid = global_hyp_pruned.look_up_table(:,i);
        index_remove = ~ismember(1:length(local_hyp_pruned{1,i}.p_exist),index_valid);
        local_hyp_pruned{1,i}.p_exist(:,index_remove) = [];
        local_hyp_pruned{1,i}.Gaussian_mixture(:,index_remove) = [];
%         local_hyp_pruned{1,i}.mean(:,index_remove) = [];
%         local_hyp_pruned{1,i}.covariance(:,:,index_remove) = [];
        local_hyp_pruned{1,i}.log_weight(:,index_remove) = [];
        
        if sum(index_remove) > 0
            for j =1:length(index_valid)
                global_hyp_pruned.look_up_table(j,i) = global_hyp_pruned.look_up_table(j,i) -sum(index_remove(1:global_hyp_pruned.look_up_table(j,i)));
            end
        end
    end
    
    % remove local hypothesis according to existence_threshold
    index_remove = [];
    revise = [];
    for i = 1:length(local_hyp_pruned)
        remove = local_hyp_pruned{1,i}.p_exist < existence_threshold;
        if sum(remove)==length(local_hyp_pruned{1,i}.p_exist)
            index_remove = [index_remove,i];
        elseif sum(remove) > 0
            revise = [revise,i];
            local_hyp_pruned{1,i}.p_exist(:,remove) = [];
            local_hyp_pruned{1,i}.Gaussian_mixture(:,remove) = [];
%             local_hyp_pruned{1,i}.mean(:,remove) = [];
%             local_hyp_pruned{1,i}.covariance(:,:,remove) = [];
            local_hyp_pruned{1,i}.log_weight(:,remove) = [];
            for j =1:size(global_hyp_pruned.look_up_table,1)
                global_hyp_pruned.look_up_table(j,i) = global_hyp_pruned.look_up_table(j,i) -sum(remove(1:global_hyp_pruned.look_up_table(j,i)));
            end
        end
    end
    local_hyp_pruned(:,index_remove) = [];
    global_hyp_pruned.look_up_table(:,index_remove) = [];
    
    if (~isempty(revise))
        [global_hyp_pruned.look_up_table, ~, dup_index] = unique(global_hyp_pruned.look_up_table,'rows');
        weight_new = zeros(1,size(global_hyp_pruned.look_up_table,1));
        for i =1:size(global_hyp_pruned.look_up_table,1)
            index_sum = (dup_index==i);
            weight_new(i) = sum(global_hyp_pruned.weight(index_sum));
        end
        global_hyp_pruned.weight=weight_new./(sum(weight_new));
    end
    
    [~,index_sorted] = sort(global_hyp_pruned.weight,'descend');
    global_hyp_pruned.weight= global_hyp_pruned.weight(index_sorted);
    global_hyp_pruned.look_up_table = global_hyp_pruned.look_up_table(index_sorted,:);
    
end