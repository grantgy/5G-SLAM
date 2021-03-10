function [result] = estimation(local_hyp,global_hyp)
    [~,index] = max(global_hyp.weight);
    data_asso = global_hyp.look_up_table(index,:);
    result = cell(size(data_asso));
    for i = 1:size(data_asso,2)
        if data_asso(1,i) == 0 || local_hyp{1,i}.p_exist(1,data_asso(1,i)) ~= 1
            result{1,i} = [];
        else
            result{1,i} = local_hyp{1,i}.mean(:,data_asso(1,i));
        end
    end
end
