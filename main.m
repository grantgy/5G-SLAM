clear all
close all

addpath(fullfile(pwd,'Functions'));

BS = [0;0;40];
VA = [200,-200,0,0;
    0,0,200,-200;
    40,40,40,40];

SP = [65,-65,-65,65;
    65,65,-65,-65];
SP=[SP;40*rand(1,4)];


N = 1;  % number of particles
K = 40;   % steps
Q = diag([0.2,0.2,0.2*pi/180,0.2].^2);  %covariance matrix of process noise
Q0= diag([0.3,0.3,0.3*pi/180,0.3].^2);    %covariance matrix to generate the particles
R = diag([0.1,0.01,0.01,0.01,0.01].^2); % covariance matrix of measurement noise
T = 0.5;  % sampling time
k_best = 20;

% initalization of vechicle 1
s1_in = [70.7285,0,pi/2,300]';  % initial state
s1_v = 22.22;                   % initial velocity
s1_ang_v = pi/10;               % initial angular velocity

% initialization of vehicle 2
s2_in = [-70.7285,0,pi/2,300]'; % initial state
s2_v = -22.22;                  % initial velocity
s2_ang_v = pi/10;               % initial angular velocity

map_BS = GM(1,[0,0,40]',eye(3));

% create the measurement
[measurement1,v1real_state]=create_measurement(s1_in,K,s1_v,s1_ang_v,T,R,BS,VA,SP);
[measurement2,v2real_state]=create_measurement(s2_in,K,s2_v,s2_ang_v,T,R,BS,VA,SP);

% GOSPA
GOSPA_p = 2;
GOSPA_c = 20;
GOSPA_alpha =2;



rmax = 200;                     % maximum sensoring range
clutter = 1/(4*rmax*pi^4);      % clutter

local_global_in = global_hypothesis(1,1); % initial global hypothesis, the base station is already know
PPP_in = [UM([],[]);UM([],[])]; % initial PPP
local_hyp_in{1,1} = MB(1,GM(1,[0,0,40]',eye(3)),0,0); %initial local hypothesis, there is the known base station

PPP_birth = [UM(0.000015,1/40000);UM(0.000015,1/40000)];

pro_sur = [0.99;0.99];
pro_det = [0.9;0.9;0.9];
% parameters for pruning
T_pruning = 0.001;       % treshold for global hypothesis
T_pruningPois=0.00001; % treshold for Possion part


Nhyp_max = 200;        % maximum number of global hypothesis
existence_threshold = 0.00001; % treshold of existence probability for Bernoulli
T_gating = 2^9;  %treshold for gating


% initialization for particles
if N >1  % if there is only one particle, that is just for testing
    s1_particle = s1_in  + mvnrnd([0;0;0;0],Q0,N)';
    s2_particle = s2_in +  mvnrnd([0;0;0;0],Q0,N)';
else
    s1_particle = s1_in ;
    s2_particle = s2_in ;
end
s1_logweight = log(ones(1,N)./N);
s2_logweight = log(ones(1,N)./N);

s1_PPP =cell(1,N);
s1_local_hyp =cell(1,N);
s1_global_hyp =cell(1,N);

s2_PPP =cell(1,N);
s2_local_hyp =cell(1,N);
s2_global_hyp =cell(1,N);

for i =1:N
    s1_PPP{1,i} = PPP_in;
    s1_local_hyp{1,i} = local_hyp_in;
    s1_global_hyp{1,i} = local_global_in;
    
    s2_PPP{1,i} = PPP_in;
    s2_local_hyp{1,i} = local_hyp_in;
    s2_global_hyp{1,i} = local_global_in;
end

v1_state = mean(s1_particle,2);
v2_state = mean(s2_particle,2);

map_v1{1,1} = MAP(map_BS,GM([],[],[]),GM([],[],[]));
map_v2{1,1} = MAP(map_BS,GM([],[],[]),GM([],[],[]));


for i = 1:K
    % prediction for particles
    s1_particle = PF_prediction(s1_particle,s1_logweight,s1_v,s1_ang_v,T,Q);
    s2_particle = PF_prediction(s2_particle,s2_logweight,s2_v,s2_ang_v,T,Q);
    
    s1_log_w_up = zeros(1,N); 
    s2_log_w_up = zeros(1,N);
    for j=1:N
        % prediction of map for each particles
        [s1_PPP{1,j},s1_local_hyp{1,j}] = PMBM_pre(s1_PPP{1,j},PPP_birth,s1_local_hyp{1,j},pro_sur);
        [s2_PPP{1,j},s2_local_hyp{1,j}] = PMBM_pre(s2_PPP{1,j},PPP_birth,s2_local_hyp{1,j},pro_sur);
        
        % update map for each particles
        [s1_PPP{1,j},s1_local_hyp{1,j},s1_global_hyp{1,j},s1_log_w_up(1,j)] = PMBM_update(s1_PPP{1,j},s1_local_hyp{1,j},clutter,pro_det, R, measurement1{1,i},s1_global_hyp{1,j},s1_particle(:,j),k_best,i);
        [s1_PPP{1,j},s1_local_hyp{1,j},s1_global_hyp{1,j}] = pruning(s1_PPP{1,j},s1_local_hyp{1,j},s1_global_hyp{1,j},T_pruning,T_pruningPois,Nhyp_max,existence_threshold);
        
        [s2_PPP{1,j},s2_local_hyp{1,j},s2_global_hyp{1,j},s2_log_w_up(1,j)] = PMBM_update(s2_PPP{1,j},s2_local_hyp{1,j},clutter,pro_det, R, measurement2{1,i},s2_global_hyp{1,j},s2_particle(:,j),k_best,i);
        [s2_PPP{1,j},s2_local_hyp{1,j},s2_global_hyp{1,j}] = pruning(s2_PPP{1,j},s2_local_hyp{1,j},s2_global_hyp{1,j},T_pruning,T_pruningPois,Nhyp_max,existence_threshold);
    end
    
    % update for particles
    [s1_particle,s1_weight,s1_state] = PF_update(s1_particle,s1_logweight,s1_log_w_up);
    [s2_particle,s2_weight,s2_state] = PF_update(s2_particle,s2_logweight,s2_log_w_up);
    
    v1_state = [v1_state,s1_state];
    v2_state = [v2_state,s2_state];
    map_v1{1,1+i} = mapreconstruction(s1_local_hyp,s1_weight,T_gating,0.1,map_BS,s1_global_hyp);
    map_v2{1,1+i} = mapreconstruction(s2_local_hyp,s2_weight,T_gating,0.1,map_BS,s2_global_hyp);
    
    [s1_particle,s1_logweight] = resampling(s1_particle,s1_state,s1_weight);
    [s2_particle,s2_logweight] = resampling(s2_particle,s2_state,s2_weight);

end

%%
T_VA=0.4;T_SP=0.2;

[error_location_v1,error_heading_v1,error_bias_v1,GOSPA_VA_v1,GOSPA_SP_v1] =Evaluation(v1_state,v1real_state,map_v1,BS,VA,SP,GOSPA_p, GOSPA_c, GOSPA_alpha,T_VA,T_SP);
[error_location_v2,error_heading_v2,error_bias_v2,GOSPA_VA_v2,GOSPA_SP_v2] =Evaluation(v2_state,v2real_state,map_v2,BS,VA,SP,GOSPA_p, GOSPA_c, GOSPA_alpha,T_VA,T_SP);

close all
figure (1)

hold on
plot(1:size(map_v1,2)-1,GOSPA_VA_v1(1,2:end),'b','LineWidth', 2)
plot(1:size(map_v1,2)-1,GOSPA_SP_v1(1,2:end),'r','LineWidth', 2)
xl=xlabel('time step');
yl=ylabel('GOSPA value' );
ll=legend('VA','SP');
title('GOSPA for Vehicle 1')
set(xl,'Interpreter','latex','FontSize',12);
set(yl,'Interpreter','latex','FontSize',12);
set(ll,'Interpreter','latex','FontSize',11,'Location','NorthEast');
pbaspect([3 1 1 ])
hold off

figure (2)

hold on
plot(1:size(map_v2,2)-1,GOSPA_VA_v2(1,2:end),'b','LineWidth', 2)
plot(1:size(map_v2,2)-1,GOSPA_SP_v2(1,2:end),'r','LineWidth', 2)
xl=xlabel('time step');
yl=ylabel('GOSPA value' );
ll=legend('VA','SP');
title('GOSPA for Vehicle 2')
set(xl,'Interpreter','latex','FontSize',12);
set(yl,'Interpreter','latex','FontSize',12);
set(ll,'Interpreter','latex','FontSize',11,'Location','NorthEast');
pbaspect([3 1 1 ])
hold off

