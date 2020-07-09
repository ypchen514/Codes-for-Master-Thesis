clear
clc
close all

sys_par = [1000,0.0001,1e4];

for s = 1:size(sys_par,1)

%% Parameters
% fixed system parameters
alpha = 1;
% tunable system parameters
p1 = -2; 
p2 = -0.5;
p3 = -2;
p4 = -0.5;

% define design variables
x1_lb = 0; 
x1_ub = 1;
x1_delta = x1_ub-x1_lb;
x2_lb = 0;
x2_ub = 1;
x2_delta = x2_ub-x2_lb;
x3_lb = 0; 
x3_ub = 1;
x3_delta = x3_ub-x3_lb;
x4_lb = 0; 
x4_ub = 1;
x4_delta = x4_ub-x4_lb;
% uncertainty
N_s = sys_par(s,1); % take N times sampling
N = sys_par(s,3);% N times sampling for kriging sensitivity analysis
tic 
%% sampling for kriging 
Q = sobolset(4,'skip',200000);% create a different sobol sequence
Q_ss_1 = net(Q,N_s);
input = Q_ss_1;% N_s by 4 matrix as input for 4 variables and N_s sampling
output = g(alpha,p1,p2,p3,p4,input)';
krig_sampling_time = toc
%% fitting
lb = [0,0,0,0];
ub = [1,1,1,1];
fittype = 1;
SCFtype = 1;
max_variance = 1;
run_num = 0;

%% Infill sampling
while max_variance > sys_par(s,2)
    run_num = run_num+1;
    % find next
    kparam = f_variogram_fit(input, output, lb, ub);
    [x_next, variance_next] = f_find_kriging_max_variance(kparam);
    output_next = g(alpha,p1,p2,p3,p4,x_next)';
    max_variance = abs(max(variance_next));
    % recurrent
    input = [input;x_next];
    output = [output;output_next];
end
total_krig_sample = length(output);
EGO_time = toc
%% Sensitivity Analysis
% fp_krig = f_predictkrige(xp, kparam);
Q = sobolset(8,'skip',20000000);% create a different sobol sequence
Q_ss = net(Q,N);
% Q_ss = lhsdesign(N,8,'smooth','off');
% load('Q_ss');

P1 = x1_lb+(Q_ss(:,1)*x1_delta);
P2 = x2_lb+(Q_ss(:,2)*x2_delta);
P3 = x3_lb+(Q_ss(:,3)*x3_delta);
P4 = x4_lb+(Q_ss(:,4)*x4_delta);
Q1 = x1_lb+(Q_ss(:,5)*x1_delta);
Q2 = x2_lb+(Q_ss(:,6)*x2_delta);
Q3 = x3_lb+(Q_ss(:,7)*x3_delta);
Q4 = x4_lb+(Q_ss(:,8)*x4_delta);

GSA_sampling_time = toc 

for j = 1:1
%% Calculation
P = [P1, P2, P3, P4];
Q = [Q1, Q2, Q3, Q4];
a = 3;

pack_size = 1e4; 
N_pack = floor(N/pack_size);
parfor i = 1:N_pack
%    if i <= N_pack
     g_P_pack = f_predictkrige(P((i-1)*pack_size+1:i*pack_size,:), kparam);
     g_Q_pack = f_predictkrige(Q((i-1)*pack_size+1:i*pack_size,:), kparam);
     g_P(i) = struct('no',i,'value',g_P_pack);
     g_Q(i) = struct('no',i,'value',g_Q_pack);
end

temp_P = [];
temp_Q = [];
for i = 1:N_pack
    temp_P = [temp_P;g_P(i).value];
    temp_Q = [temp_Q;g_Q(i).value];
end
g_P = temp_P';
g_Q = temp_Q';

sub_time = toc
var_g = 1/N*(sum(g_P.^2))-(1/N*sum(g_P))^2; %Total variance 

for i = 1:4
    PP = P;
    QQ = Q;
    PP(:,i) = QQ(:,i);
    R(i) = struct('number',i,'in',PP); % will generate R matrix as sample matrix
    g_R(i) = struct('no',i,'value',g(alpha,p1,p2,p3,p4,R(i).in));% calulate each g(R^j)
    C(i) = struct('no',i,'value',1/N*(sum(g_Q.*(g_R(i).value-g_P))));
    D(i) = struct('no',i,'value',1/(2*N)*sum((g_P-g_R(i).value).^2));
    Sj(i) = struct('design_variable',i,'main_effect',C(i).value/var_g,'total_effect',D(i).value/var_g);
    
end    

Sj_s(j) = struct('no',j,'S1',Sj(1).main_effect,'S2',Sj(2).main_effect,'S3',Sj(3).main_effect,'S4',Sj(4).main_effect);

end
total_time = toc

% simulation data
Sim_data(s) = struct('Ini_sp',sys_par(s,1),'max_var',sys_par(s,2),'MCS_sp',sys_par(s,3),'Krig_samp_t',krig_sampling_time,'EGO_N', total_krig_sample,'EGO_t',EGO_time,'GSA_sam_t',GSA_sampling_time,'total_t',total_time,'Sj',Sj)

clear g_P;
clear g_Q;
end






