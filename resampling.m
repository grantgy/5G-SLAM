function [particle,log_w_new] = resampling(particle,State,weight)
N = size(particle,2);
if N>1
    particle = State + mvnrnd([0;0;0;0],diag([0.3,0.3,0.3*pi/180,0.3].^2),N)';
    log_w_new = log(ones(1,N)./N);
else
    log_w_new = log(weight);

end