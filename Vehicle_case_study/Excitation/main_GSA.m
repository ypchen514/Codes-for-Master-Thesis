%% SOBOL INDICES ON GLOBAL SENSITIVITY FOR TADPOLE THREE WHEEL VEHICLE
clc
close all

clear
total_count = 0;
% for iii = 0:2:8
%     for jjj = 0:2:8
total_count = total_count+1;
%% 1. Complex Model Sampling
% In this sampling stage, we first need to outcome the input and output
% data as the learing data for kriging fitting.
% input x: x is derived from sobol sequence depends on variable numbers and
% sampling numbers.
% output y: is the total squared error of output track and nominal values
% set.

runcycle = 8; % total sampling times
n_data_collect = 10000; % cut the output track into # pieces
tic
%---Nominal Set------------------------------------------------------
% give maneuver parameters
case_name = 'circle';
% sys_par = [2,2,0,1,1]; %need a1, a2, b(better = 0), zoom_rate, zoom_length
sys_par = 5;
sim_time = 40;
whole_par; % generate the nominal value of system parameters
% Drive once to generate maneuver command
[waypoint, require_velocity] = waypoints(sys_par,case_name); % create path following points
sim('tadpole_dynamic_close_loop'); % drive once, ouput its track's x and y points in cartisian coordinate
x_nominal = interp(linspace(0,length(x.Data),length(x.Data)),x.Data,n_data_collect);% x of each points
y_nominal = interp(linspace(0,length(y.Data),length(y.Data)),y.Data,n_data_collect);% y of each points
% Drive command
steer_command = steer_com.Data(1,:);
drive_command = drive_com.Data(1,:);

%---Complex System Sampling------------------------------------------
% Create parameter sets from low discrepency sampling
[I_x,I_y,I_z,C_alpha,C_beta,SAP,R_0,k_z,delta_fL,delta_fR,V_wind,rho,Cd,FA,mu_rolling_1,mu_rolling_2,mu_dp,mu_lp,Vx,Vy,Vz,DOE] = rand_par(runcycle);
% insert parameters into simulink model
% take the normalized number into design parameters of vehicle
fraction = 4;% each fraction has 100 samples
total_sq_err = [];
x_sq_err = [];
y_sq_err = [];
for j = 1:(runcycle/fraction)
    for i = (1+fraction*(j-1)):fraction*j
        in(i) = Simulink.SimulationInput('tadpole_dynamic_open_loop1');%drive open loop vehicle
        in(i) = in(i).setVariable('l_1',l_1);
        in(i) = in(i).setVariable('l_2',l_2);
        in(i) = in(i).setVariable('w_1',w_1);
        in(i) = in(i).setVariable('w_2',w_2);
        in(i) = in(i).setVariable('h',h);
        in(i) = in(i).setVariable('p',p);
        in(i) = in(i).setVariable('m',m);
        in(i) = in(i).setVariable('g',9.81);
        in(i) = in(i).setVariable('tire_r',0.3302);
        in(i) = in(i).setVariable('tire_wdt',0.025);
        in(i) = in(i).setVariable('I_x',I_x(i));
        in(i) = in(i).setVariable('I_y',I_y(i));
        in(i) = in(i).setVariable('I_z',I_z(i));
        in(i) = in(i).setVariable('C_alpha',C_alpha(i));
        in(i) = in(i).setVariable('C_beta',C_beta(i));
        in(i) = in(i).setVariable('SAP',SAP(i));
        in(i) = in(i).setVariable('R_0',R_0(i));
        in(i) = in(i).setVariable('k_z',k_z(i));
        in(i) = in(i).setVariable('delta_fL',delta_fL(i));
        in(i) = in(i).setVariable('delta_fR',delta_fR(i));
        in(i) = in(i).setVariable('V_wind',V_wind(i));
        in(i) = in(i).setVariable('rho',rho(i));
        in(i) = in(i).setVariable('Cd',Cd(i));
        in(i) = in(i).setVariable('FA',FA(i));
        in(i) = in(i).setVariable('mu_rolling_1',mu_rolling_1(i));
        in(i) = in(i).setVariable('mu_rolling_2',mu_rolling_2(i));
        in(i) = in(i).setVariable('mu_dp',mu_dp(i));
        in(i) = in(i).setVariable('mu_lp',mu_lp(i));
        in(i) = in(i).setVariable('Vx',0);
        in(i) = in(i).setVariable('Vy',0);
        in(i) = in(i).setVariable('Vz',0);
        in(i) = in(i).setVariable('I_f_wheel',0.3);
        in(i) = in(i).setVariable('I_r_wheel',0.5);
        in(i) = in(i).setVariable('drive_command',drive_command);
        in(i) = in(i).setVariable('steer_command',steer_command);
        in(i) = in(i).setVariable('sim_time',sim_time);
    end
    % set_param('tadpole_dynamic_open_loop1','AlgebraicLoopSolver','TrustRegion')
    out = parsim(in,'ShowSimulationManager','on','ShowProgress','on');%parellel computation
    par_toc = toc;
    
    for i = (1+fraction*(j-1)):fraction*j
        x_open_pt(i,:) = interp(linspace(0,length(out(1,i).x_open.Data),length(out(1,i).x_open.Data)),out(1,i).x_open.Data,n_data_collect);
        y_open_pt(i,:) = interp(linspace(0,length(out(1,i).y_open.Data),length(out(1,i).y_open.Data)),out(1,i).y_open.Data,n_data_collect);
        plot(x_open_pt(i,:),y_open_pt(i,:));
        hold on
        %calculate sum error
        [total_sq_err_temp,x_sq_err_temp,y_sq_err_temp] = sq_error(x_open_pt,y_open_pt,x_nominal,y_nominal,n_data_collect); %calculate square error
    end
    total_sq_err = [total_sq_err;total_sq_err_temp];
    x_sq_err = [x_sq_err;x_sq_err_temp];
    y_sq_err = [y_sq_err;y_sq_err_temp];
end
interp_toc = toc;
%% 2. Kriging Fitting and Infill Sampling
% Initial fitting setting
input = DOE*2-1;
output = total_sq_err;
lb = -ones(1,size(input,2));
ub = ones(1,size(input,2));
fittype = 1;
SCFtype = 1;
max_variance = 1;
run_num = 0;
% not doing infill
kparam = f_variogram_fit(input, output, lb, ub);% Fitting of kriging with i/o
[x_next, variance_next] = f_find_kriging_max_variance(kparam);
% Infill sampling
% while max_variance > 0.0008 && run_num < 500
%     run_num = run_num+1;% total infill sampling numbers
%     kparam = f_variogram_fit(input, output, lb, ub);% Fitting of kriging with i/o
%     [x_next, variance_next] = f_find_kriging_max_variance(kparam);
%     % use x_next to outcome next output
%     [I_x,I_y,I_z,C_alpha,C_beta,SAP,R_0,k_z,delta_fL,delta_fR,V_wind,rho,Cd,FA,mu_rolling_1,mu_rolling_2,mu_dp,mu_lp,Vx,Vy,Vz] = infill_par(x_next);
%     sim('tadpole_dynamic_open_loop',20);
%     % linear interpolation of new outcome
%     x_open_EGO_pt = interp(linspace(0,length(x_open.Data),length(x_open.Data)),x_open.Data,n_data_collect);
%     y_open_EGO_pt = interp(linspace(0,length(y_open.Data),length(y_open.Data)),y_open.Data,n_data_collect);
%     [total_sq_err_EGO,x_sq_err_EGO,y_sq_err_EGO] = sq_error(x_open_EGO_pt,y_open_EGO_pt,x_nominal,y_nominal,n_data_collect); %calculate square error
% save data
%     plot(x_open_EGO_pt,y_open_EGO_pt);
%     x_open_pt = [x_open_pt;x_open_EGO_pt];
%     y_open_pt = [y_open_pt;y_open_EGO_pt];
%     x_sq_err = [x_sq_err;x_sq_err_EGO];
%     y_sq_err = [y_sq_err;y_sq_err_EGO];
%     total_sq_err = [total_sq_err;total_sq_err_EGO];
%     input = [input;x_next];
%     output = total_sq_err;
%     max_variance = abs(variance_next);
%     variance_data(run_num) = variance_next;
% end
% EGO_toc = toc;


%% 3. Sensitivity Analysis

N = 5e5; % GSA monte carlo sampling numbers
Q = sobolset(12,'skip',20000000);% create a different sobol sequence
Q_ss = net(Q,N);
P1 = (Q_ss(:,1)*2-1);
P2 = (Q_ss(:,2)*2-1);
P3 = (Q_ss(:,3)*2-1);
P4 = (Q_ss(:,4)*2-1);
P5 = (Q_ss(:,5)*2-1);
P6 = (Q_ss(:,6)*2-1);
Q1 = (Q_ss(:,7)*2-1);
Q2 = (Q_ss(:,8)*2-1);
Q3 = (Q_ss(:,9)*2-1);
Q4 = (Q_ss(:,10)*2-1);
Q5 = (Q_ss(:,11)*2-1);
Q6 = (Q_ss(:,12)*2-1);

GSA_sampling_time = toc

% Calculation
P = [P1, P2, P3, P4, P5, P6];
Q = [Q1, Q2, Q3, Q4, Q5, Q6];
g_P = [];
g_Q = [];
for k = 1:length(P)/10000
    g_P_temp = f_predictkrige(P((k-1)*1000+1:k*1000,:), kparam)';
    g_Q_temp= f_predictkrige(Q((k-1)*1000+1:k*1000,:), kparam)';
    g_P = [g_P g_P_temp];
    g_Q = [g_Q g_Q_temp];
end

sub_time = toc
var_g = 1/N*(sum(g_P.^2))-(1/N*sum(g_P))^2; %Total variance


for i = 1:6
    PP = P;
    QQ = Q;
    PP(:,i) = QQ(:,i);
    R(i) = struct('number',i,'in',PP); % will generate R matrix as sample matrix
    g_RR = [];
    for k = 1:length(PP)/10000
        g_RR_temp = f_predictkrige(PP((k-1)*1000+1:k*1000,:), kparam)';
        g_RR = [g_RR g_RR_temp];
    end
    g_R(i) = struct('no',i,'value',g_RR);% calulate each g(R^j)
    C(i) = struct('no',i,'value',1/N*(sum(g_Q.*(g_R(i).value-g_P))));
    D(i) = struct('no',i,'value',1/(2*N)*sum((g_P-g_R(i).value).^2));
    Sj(i) = struct('design_variable',i,'main_effect',C(i).value/var_g,'total_effect',D(i).value/var_g);
    
end
MCS_time = toc
%     Sj_s(total_count) = struct('no',total_count,'i',i,'j',j,'MCS_time',MCS_time,'Sj',Sj);
%
%     end
% end

