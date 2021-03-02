function [PPP,local_hyp] = PMBM_pre(PPP,PPP_birth,local_hyp,pro_sur)
    % predicted PPP intensity for surriving undetected objects
    %PPP.weight = PPP.weight.*pro_sur;
    for i =1 : length(PPP)
        PPP(i,1).weight = PPP(i,1).weight.*pro_sur(i,1);
        PPP(i,1) = UM([PPP(i,1).weight,PPP_birth(i,1).weight],[PPP(i,1).pro,PPP_birth(i,1).pro]);
    end 
    
    if length(local_hyp) > 1
        for i = 2:length(local_hyp)
            for j =1:length(local_hyp{1,i}.p_exist)
                local_hyp{1,i}.p_exist(j) = local_hyp{1,i}.p_exist(j).*(sum(local_hyp{1,i}.Gaussian_mixture(j).weight.*pro_sur(1:size(local_hyp{1,i}.Gaussian_mixture(j).weight,2),1)'));
            end
        end
    end
%     
%     for i = 1: length(PPP.weight)
%         [PPP.mean(:,i),PPP.covariance(:,:,i)] = kf_prediction(PPP.mean(:,i),PPP.covariance(:,:,i),F,Q);
%     end
%     
    % predicted PPP intensity for undetected objects (the surrival and birth)
%     PPP = GM([PPP.weight,PPP_birth.weight],[PPP.mean,PPP_birth.mean],cat(3,PPP.covariance,PPP_birth.covariance));
    
    % predicted MBM for detected objects
    
%     for i = 1:length(local_hyp)
%        local_hyp{1,i}.p_exist = local_hyp{1,i}.p_exist.*sum(pro_sur);
%        for j = 1:length(local_hyp{1,i}.p_exist)
%            [local_hyp{1,i}.mean(:,j),local_hyp{1,i}.covariance(:,:,j)] = kf_prediction(local_hyp{1,i}.mean(:,j),local_hyp{1,i}.covariance(:,:,j),F,Q);
%        end
%     end    
end