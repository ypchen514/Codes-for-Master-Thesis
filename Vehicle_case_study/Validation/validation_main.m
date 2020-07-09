%% SOBOL INDICES ON GLOBAL SENSITIVITY FOR TADPOLE THREE WHEEL VEHICLE -Validation
clc
close all

clear
tic
total_count = 0;
sys_par_data = [7,0,0,1,1];

sys_par = sys_par_data;
iii = sys_par(1);
jjj = sys_par(2);

total_count = total_count+1;
%% 1. Complex Model Sampling
% In this sampling stage, we first need to outcome the input and output
% data as the learing data for kriging fitting.
% input x: x is derived from sobol sequence depends on variable numbers and
% sampling numbers.
% output y: is the total squared error of output track and nominal values
% set.

runcycle = 500; % total sampling times
n_data_collect = 100000; % cut the output track into # pieces
tic
%---Nominal Set------------------------------------------------------
% give maneuver parameters
sys_par = [iii,jjj,0,1,1]; %need a1, a2, b(better = 0), zoom_rate, zoom_length
%         sys_par = [iii,jjj];
sim_time = 85;
whole_par; % generate the nominal value of system parameters
% Drive once to generate maneuver command
%     [waypoint, require_velocity] = waypoints(sys_par,case_name); % create path following points
[waypoint, require_velocity] = validation_waypoint(sys_par);
sim('tadpole_dynamic_close_loop'); % drive once, ouput its track's x and y points in cartisian coordinate
x_nominal = interp(linspace(0,length(x.Data),length(x.Data)),x.Data,n_data_collect);% x of each points
y_nominal = interp(linspace(0,length(y.Data),length(y.Data)),y.Data,n_data_collect);% y of each points
% Drive command
steer_command = steer_com.Data(1,:);
drive_command = drive_com.Data(1,:);
plot(x_nominal,y_nominal);

%---Complex System Sampling------------------------------------------
% Create parameter sets from low discrepency sampling
[I_x,I_y,I_z,C_alpha,C_beta,SAP,R_0,k_z,delta_fL,delta_fR,V_wind,rho,Cd,FA,mu_rolling_1,mu_rolling_2,mu_dp,mu_lp,Vx,Vy,Vz,DOE] = rand_par3(runcycle);
% insert parameters into simulink model
% take the normalized number into design parameters of vehicle
fraction = runcycle;% each fraction has 200 samples
j=1;

for i = 1:fraction
    in(i) = Simulink.SimulationInput('tadpole_dynamic_open_loop2');%drive open loop vehicle
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
    in(i) = in(i).setVariable('I_x',I_x(i+(j-1)*fraction));
    in(i) = in(i).setVariable('I_y',I_y(i+(j-1)*fraction));
    in(i) = in(i).setVariable('I_z',I_z(i+(j-1)*fraction));
    in(i) = in(i).setVariable('C_alpha',C_alpha(i+(j-1)*fraction));
    in(i) = in(i).setVariable('C_beta',C_beta(i+(j-1)*fraction));
    in(i) = in(i).setVariable('SAP',SAP(i+(j-1)*fraction));
    in(i) = in(i).setVariable('R_0',R_0(i+(j-1)*fraction));
    in(i) = in(i).setVariable('k_z',k_z(i+(j-1)*fraction));
    in(i) = in(i).setVariable('delta_fL',delta_fL(i+(j-1)*fraction));
    in(i) = in(i).setVariable('delta_fR',delta_fR(i+(j-1)*fraction));
    in(i) = in(i).setVariable('V_wind',V_wind(i+(j-1)*fraction));
    in(i) = in(i).setVariable('rho',rho(i+(j-1)*fraction));
    in(i) = in(i).setVariable('Cd',Cd(i+(j-1)*fraction));
    in(i) = in(i).setVariable('FA',FA(i+(j-1)*fraction));
    in(i) = in(i).setVariable('mu_rolling_1',mu_rolling_1(i+(j-1)*fraction));
    in(i) = in(i).setVariable('mu_rolling_2',mu_rolling_2(i+(j-1)*fraction));
    in(i) = in(i).setVariable('mu_dp',mu_dp(i+(j-1)*fraction));
    in(i) = in(i).setVariable('mu_lp',mu_lp(i+(j-1)*fraction));
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
for i = 1:fraction
    % Real Output Data
    x_open_pt(i,:) = interp(linspace(0,length(out(i).x_open.signals.values),length(out(i).x_open.signals.values)),out(i).x_open.signals.values,n_data_collect);
    y_open_pt(i,:) = interp(linspace(0,length(out(i).y_open.signals.values),length(out(i).y_open.signals.values)),out(i).y_open.signals.values,n_data_collect);
    ax_open_pt(i,:) = interp(linspace(0,length(out(i).a_x.signals.values),length(out(i).a_x.signals.values)),out(i).a_x.signals.values,n_data_collect);
    ay_open_pt(i,:) = interp(linspace(0,length(out(i).a_y.signals.values),length(out(i).a_y.signals.values)),out(i).a_y.signals.values,n_data_collect);
    yaw_open_pt(i,:) = interp(linspace(0,length(out(i).yaw_rate.signals.values),length(out(i).yaw_rate.signals.values)),out(i).yaw_rate.signals.values,n_data_collect);
    time_pt(i,:) = interp(linspace(0,length(out(i).x_open1.time),length(out(i).x_open1.time)),out(i).x_open1.time,n_data_collect);
    plot(x_open_pt(i,:),y_open_pt(i,:));
    hold on
end
