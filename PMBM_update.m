function [PPP_update,local_hyp_update,global_hyp_update,log_w_sum] = PMBM_update(PPP,local_hyp,clutter,pro_det, R, measurement,globH,particle,k_best,time)


    N_object = size(local_hyp,2);
    local_hyp_update = [local_hyp,cell(1,size(measurement,2))];
    % undetected objects that remain undetected (PPP)
    PPP_update = PPP;
    for i=1:length(PPP)
        PPP_update(i,1).weight = PPP(i,1).weight.* (1-pro_det(i+1,1));
    end
    
    % object detected for the first time
    
    %logweig_newB = zeros(1,size(measurement,2));

    
    for j = 1:size(measurement,2)
        ro = zeros(length(PPP),1);
        GM_new1 = cell(length(PPP),1);
        for m = 1:length(PPP)
            GM_new = GM([],[],[]);% PPP(m,1);
            for i = 1:length(PPP(m,1).weight)
                [mean, covariance] = newBern(measurement(:,j),9*R,particle,m+1);
                
                [~, ~, ~,R_inv] = CKF(particle,measurement(:,j),mean, covariance, 9*R, m+1);

                %local_hyp_update{1,k}.his = [local_hyp_update{1,k}.his,j];
                z_inv = measurement_function(mean,particle,[0,0,40]',m+1,measurement(:,j));                      
                nor_factor = pro_det(m+1,1)*mvnpdf(measurement(:,j),z_inv,R_inv)*PPP(m,1).weight(i);
                %nor_factor = pro_det(m+1,1);
                
                
                GM_new.mean = [GM_new.mean,mean];
                GM_new.covariance = cat(3,GM_new.covariance,covariance);
                
                
%               [GM_new.mean(:,i), GM_new.covariance(:,:,i), z_inv, R_inv] = CKF(particle,measurement(:,j), PPP(m,1).mean(:,i), PPP(m,1).covariance(:,:,i), R, m);
%                nor_factor = pro_det(m,1)*PPP(m,1).weight(i)*mvnpdf(measurement(:,j),z_inv, R_inv);
                %GM_new.weight(i) = nor_factor;
                GM_new.weight = [GM_new.weight,nor_factor];
                ro(m,1) = ro(m,1)+ nor_factor;
            end
            GM_new1{m,1} = GM_new;
        end
        %GM_new.weight = GM_new.weight ./ ro;
        %[x_post,P_post] = merge(GM_new); 
        %ro = sum(ro);
        p_new = sum(ro) / (sum(ro)+clutter)/10;
%         x_post1 = cell(length(PPP),1);
%         P_post1 = cell(length(PPP),1);
        x_post1 = [];
        P_post1 = [];
        for m =1:length(PPP)
            GM_new1{m,1}.weight = GM_new1{m,1}.weight ./ ro(m,1);
            [x_post,P_post] = merge(GM_new1{m,1});
            x_post1 = [x_post1,x_post];
            P_post1 = cat(3,P_post1,P_post);
        end
        
        %logweig_newB(1,j) = log(ro+clutter);
        local_hyp_update{1,N_object + j} = MB(p_new,GM(ro'./sum(ro),x_post1,P_post1),log(sum(ro)+clutter),time);
        %local_hyp_update{1,N_object + j}.time = t;   need to be refined
        %local_hyp_update{1,N_object + j}.his = j;
        
    end
    

    
    % update the local hypothesis(leaf) for each MB
    
    % update the base station
    local_hyp_update{1,1} = MB([],[],[],0);
    for i =1:length(local_hyp{1,1}.p_exist)
        
        pro_d = pro_det(1,1);
        
        % misdetected
        pro = 0;
        New_GM = GM([],[],[]);
        New_GM.weight = 1-pro_d;
        New_GM.mean = local_hyp{1,1}.Gaussian_mixture(i).mean(:,1);
        New_GM.covariance = local_hyp{1,1}.Gaussian_mixture(i).covariance(:,:,1);
        %local_hyp_update{1,k}.his = [local_hyp_update{1,k}.his,0];
        pro = pro + (1-pro_d);
        
        New_GM.weight = New_GM.weight./pro;
        local_hyp_update{1,1}.p_exist = [local_hyp_update{1,1}.p_exist, local_hyp{1,1}.p_exist(i)*pro/(1-local_hyp{1,1}.p_exist(i)+local_hyp{1,1}.p_exist(i)*pro)];
        local_hyp_update{1,1}.Gaussian_mixture = [local_hyp_update{1,1}.Gaussian_mixture, New_GM];
        local_hyp_update{1,1}.log_weight = [local_hyp_update{1,1}.log_weight, log(1-local_hyp{1,1}.p_exist(i)+local_hyp{1,1}.p_exist(i)*pro )];

        
        % detected
        for j = 1:size(measurement,2)
            
            pro = 0;
            New_GM = GM([],[],[]);
            
            [x_post, P_poster, z_inv,R_inv] = CKF(particle,measurement(:,j),local_hyp{1,1}.Gaussian_mixture(i).mean(:,1), local_hyp{1,1}.Gaussian_mixture(i).covariance(:,:,1), R, 1);
            New_GM.weight = 1;
            New_GM.mean = [New_GM.mean, x_post];
            New_GM.covariance = cat(3,New_GM.covariance,P_poster);
            %local_hyp_update{1,k}.his = [local_hyp_update{1,k}.his,j];
            %z_inv = measurement_function(x_post,particle,[0,0,40]',1,measurement(:,j));
            %z_inv = Cali(z_inv, measurement(:,j));
            pro = pro + pro_d*mvnpdf(measurement(:,j),z_inv,16*R_inv);
            
            %New_GM.weight = New_GM.weight./pro;
            local_hyp_update{1,1}.p_exist = [local_hyp_update{1,1}.p_exist, 1];
            local_hyp_update{1,1}.Gaussian_mixture = [local_hyp_update{1,1}.Gaussian_mixture, New_GM];
            local_hyp_update{1,1}.log_weight = [local_hyp_update{1,1}.log_weight, log(local_hyp{1,1}.p_exist(i)*pro )];
            
        end
        
    end
    
    
    % uppdate VA/SP
    for k = 2: N_object
        local_hyp_update{1,k} = MB([],[],[],local_hyp{1,k}.birth_time);
        for i =1:length(local_hyp{1,k}.p_exist)
            if size(local_hyp{1,k}.Gaussian_mixture(i).mean,2) == 2
                if norm([particle(1:2,1);0]-local_hyp{1,k}.Gaussian_mixture(i).mean(:,2))<=50+2*sqrt(diag(local_hyp{1,k}.Gaussian_mixture(i).covariance(:,:,2)))
                    pro_d = pro_det(2:3,1);
                else
                    pro_d = [pro_det(2,1);0];
                end
            else
                pro_d = pro_det(2,1);
            end
            

            % misdetected
            pro = 0;
            New_GM = GM([],[],[]);
            for m = 1:length(pro_d)
                New_GM.weight = [New_GM.weight,(1-pro_d(m,1))*local_hyp{1,k}.Gaussian_mixture(i).weight(m)];
                New_GM.mean = [New_GM.mean,local_hyp{1,k}.Gaussian_mixture(i).mean(:,m)];
                New_GM.covariance = cat(3,New_GM.covariance,local_hyp{1,k}.Gaussian_mixture(i).covariance(:,:,m));
                pro = pro + (1-pro_d(m,1))*local_hyp{1,k}.Gaussian_mixture(i).weight(m);
            end
            New_GM.weight = New_GM.weight./pro;
            local_hyp_update{1,k}.p_exist = [local_hyp_update{1,k}.p_exist, local_hyp{1,k}.p_exist(i)*pro/(1-local_hyp{1,k}.p_exist(i)+local_hyp{1,k}.p_exist(i)*pro)];
            local_hyp_update{1,k}.Gaussian_mixture = [local_hyp_update{1,k}.Gaussian_mixture, New_GM];
            local_hyp_update{1,k}.log_weight = [local_hyp_update{1,k}.log_weight, log(1-local_hyp{1,k}.p_exist(i)+local_hyp{1,k}.p_exist(i)*pro )];
            

            % detected
            
            for j = 1:size(measurement,2)
                pro = 0;
                New_GM = GM([],[],[]);
                
                for m =1:length(pro_d)
                    if pro_d(m) ~=0
                        
                        [x_post, P_poster, z_inv,R_inv] = CKF(particle,measurement(:,j),local_hyp{1,k}.Gaussian_mixture(i).mean(:,m), local_hyp{1,k}.Gaussian_mixture(i).covariance(:,:,m), 9*R, m+1);

                        
                        New_GM.weight = [New_GM.weight,pro_d(m,1)*mvnpdf(measurement(:,j),z_inv,R_inv)*local_hyp{1,k}.Gaussian_mixture(i).weight(m)];
                        New_GM.mean = [New_GM.mean, x_post];
                        New_GM.covariance = cat(3,New_GM.covariance,P_poster);
                        pro = pro + New_GM.weight(m);
                    end
                end
                New_GM.weight = New_GM.weight./pro;
                local_hyp_update{1,k}.p_exist = [local_hyp_update{1,k}.p_exist, 1];
                local_hyp_update{1,k}.Gaussian_mixture = [local_hyp_update{1,k}.Gaussian_mixture, New_GM];
                local_hyp_update{1,k}.log_weight = [local_hyp_update{1,k}.log_weight, log(local_hyp{1,k}.p_exist(i)*pro )];
                
            end

        end
    end
    
      
    globH_pre =globH.look_up_table;
    globH_pre_weight=globH.weight;
    globH_new = [];
    globH_new_weight = [];
    N_mea = size(measurement,2);
    
    
    cost2 = inf*ones(size(measurement,2));
    for i =1:size(measurement,2)
        cost2(i,i) = -local_hyp_update{1,N_object+i}.log_weight(1);
    end
    
    for j = 1:length(globH_pre_weight)
        cost1 = inf*ones(size(measurement,2),N_object);
        for i =1:size(globH_pre,2)
            if globH_pre(j,i)~=0
                cost1(:,i) = (local_hyp_update{1,i}.log_weight(1+(globH_pre(j,i)-1)*(N_mea+1))-local_hyp_update{1,i}.log_weight(2+(globH_pre(j,i)-1)*(N_mea+1):globH_pre(j,i)*(N_mea+1)))';
            end
        end
        cost_matrix = [cost1,cost2];

        
        [assignments,weights]= murty(cost_matrix,k_best);
        
        log_w = log(globH_pre_weight(1,j))*ones(1,size(assignments,1));
        globH = zeros(size(assignments,1),N_object+N_mea);
        
        for i = N_object +1 : N_object+N_mea
            for k = 1:size(assignments,1)
                index_object =  find(assignments(k,:)==i,1);
                if(isempty(index_object))
                    index= 0;
                else 
                    index= 1;
                    log_w(1,k) = log_w(1,k) + local_hyp_update{1,i}.log_weight(1);
                    
                end
                globH(k,i)=index;
            end
        end
        
        for i = 1:N_object
            index_pre = globH_pre(j,i);
            for k = 1:size(assignments,1)
                index_object =  find(assignments(k,:)==i);
                
                if(isempty(index_object))
                    if index_pre ==0
                        index =0;
                    else
                        index= (index_pre-1) * (N_mea+1) +1;
                        log_w(1,k) = log_w(1,k) + local_hyp_update{1,i}.log_weight(index);
                    end
                else
                    
                    index= (index_pre-1) * (N_mea+1) +1+index_object;
                    log_w(1,k) = log_w(1,k) + local_hyp_update{1,i}.log_weight(index);
                end
                globH(k,i)=index;
            end
        end
        globH_new = [globH_new;globH];
        globH_new_weight = [globH_new_weight,log_w];

    end
        globH_new_weight_max = max(globH_new_weight);
        globH_new_weight = exp(globH_new_weight-max(globH_new_weight));
        log_w_sum = log(sum(globH_new_weight))+globH_new_weight_max;
        globH_new_weight = globH_new_weight./sum(globH_new_weight);
        global_hyp_update = global_hypothesis(globH_new_weight,globH_new);    
    
    
    
end
        