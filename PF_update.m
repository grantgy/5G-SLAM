function [particle,weight,State] = PF_update(particle,logweight,log_w)
    log_w_new = logweight + log_w;
    log_w_new = exp(log_w_new-max(log_w_new));
    weight = log_w_new./sum(log_w_new);
    
    State = sum(particle.*weight,2);
    State(3) = mod(State(3),2*pi);
%     N = size(particle,2);
%     if N>1
%         particle = State + mvnrnd([0;0;0;0],diag([0.3,0.3,0.3*pi/180,0.3].^2),N)';
%         log_w_new = log(ones(1,N)./N);
%     else
%         log_w_new = log(log_w_new);
%     end
end