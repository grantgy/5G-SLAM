function [particle,weight,State] = PF_update(particle,logweight,log_w)
    log_w_new = logweight + log_w;
    log_w_new = exp(log_w_new-max(log_w_new));
    weight = log_w_new./sum(log_w_new);
    
    State = sum(particle.*weight,2);
    State(3) = mod(State(3),2*pi);
end
