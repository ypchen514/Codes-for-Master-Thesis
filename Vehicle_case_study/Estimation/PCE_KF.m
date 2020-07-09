%% Polynomial Chaos Expansion Based Kalman Filter on Parameter Estimation of Complex system
%====================================================================================
clc
clear
close all

%% Univariate Polynomials
% Create Legendre Function
% Variable 1
syms x_1
phi_1 = [1;x_1;1/2*(3*x_1^2-1);1/2*(5*x_1^3-3*x_1);...% 0-3
    1/8*(35*x_1^4-30*x_1^2+3);1/8*(63*x_1^5-70*x_1^3+15*x_1);... % 4-5
    1/16*(231*x_1^6-315*x_1^4+105*x_1^2-5);1/16*(429*x_1^7-693*x_1^5+315*x_1^3-35*x_1)]; % 6-7
subs(phi_1,0.5);

% Variable 2
syms x_2
phi_2 = [1;x_2;1/2*(3*x_2^2-1);1/2*(5*x_2^3-3*x_2);...% 0-3
    1/8*(35*x_2^4-30*x_2^2+3);1/8*(63*x_2^5-70*x_2^3+15*x_2);... % 4-5
    1/16*(231*x_2^6-315*x_2^4+105*x_2^2-5);1/16*(429*x_2^7-693*x_2^5+315*x_2^3-35*x_2)]; % 6-7

% Variable 3
syms x_3
syms x_3
phi_3 = [1;x_3;1/2*(3*x_3^2-1);1/2*(5*x_3^3-3*x_3);...% 0-3
    1/8*(35*x_3^4-30*x_3^2+3);1/8*(63*x_3^5-70*x_3^3+15*x_3);... % 4-5
    1/16*(231*x_3^6-315*x_3^4+105*x_3^2-5);1/16*(429*x_3^7-693*x_3^5+315*x_3^3-35*x_3)]; % 6-7

% Variable 4
syms x_4
syms x_4
phi_4 = [1;x_4;1/2*(3*x_4^2-1);1/2*(5*x_4^3-3*x_4);...% 0-3
    1/8*(35*x_4^4-30*x_4^2+3);1/8*(63*x_4^5-70*x_4^3+15*x_4);... % 4-5
    1/16*(231*x_4^6-315*x_4^4+105*x_4^2-5);1/16*(429*x_4^7-693*x_4^5+315*x_4^3-35*x_4)]; % 6-7

% Variable 5
syms x_5
syms x_5
phi_5 = [1;x_5;1/2*(3*x_5^2-1);1/2*(5*x_5^3-3*x_5);...% 0-3
    1/8*(35*x_5^4-30*x_5^2+3);1/8*(63*x_5^5-70*x_5^3+15*x_5);... % 4-5
    1/16*(231*x_5^6-315*x_5^4+105*x_5^2-5);1/16*(429*x_5^7-693*x_5^5+315*x_5^3-35*x_5)]; % 6-7

% Variable 6
syms x_6
syms x_6
phi_6 = [1;x_6;1/2*(3*x_6^2-1);1/2*(5*x_6^3-3*x_6);...% 0-3
    1/8*(35*x_6^4-30*x_6^2+3);1/8*(63*x_6^5-70*x_6^3+15*x_6);... % 4-5
    1/16*(231*x_6^6-315*x_6^4+105*x_6^2-5);1/16*(429*x_6^7-693*x_6^5+315*x_6^3-35*x_6)]; % 6-7


%% Multivariate Polynomials
% The polynomial chaos expansion term are the same for input and system
% output. So, we only need to do this process once to generate

% Multiply Term
P = 5; % truncation, P is the highest order of the polynomial
phi_term1 = [0;1;2;3;4;5;6;7];
phi_term2 = [0;1;2;3;4;5;6;7];
phi_term3 = [0;1;2;3;4;5;6;7];
phi_term4 = [0;1;2;3;4;5;6;7];
phi_term5 = [0;1;2;3;4;5;6;7];
phi_term6 = [0;1;2;3;4;5;6;7];

tic
PHI = [];
location = [];
Pos = 0;
for i = 1:length(phi_term1)
    for j = 1:length(phi_term2)
        for l = 1:length(phi_term3)
            for m = 1:length(phi_term4)
                for n = 1:length(phi_term5)
                    for o = 1:length(phi_term6)
                        if phi_term1(i)+phi_term2(j)+phi_term3(l)+phi_term4(m)+phi_term5(n)+phi_term6(o) <= P
                            Mul = [i-1,j-1,l-1,m-1,n-1,o-1];
                            Pos = Pos+1;
                            order = sum(Mul);
                            if nnz(Mul) ==1
                                location = [location;[Pos,order]];
                            end
                            PHI_temp = phi_1(i)*phi_2(j)*phi_3(l)*phi_4(m)*phi_5(n)*phi_6(o);
                            PHI = [PHI;PHI_temp]; % final multivariate polynomial equation
                        end
                    end
                end
            end
        end
    end
end
toc

%% Collocation Point
Qs = sobolset(6,'skip',30);
Q = net(Qs,3*length(PHI)); % 0~1

%% Coefficient
% By giving random variable(-1~1) with n dimension, we can use polynomial
% chaos expansion to infer what system input is.
% Define x
x = (Q-0.5)*2; % -1 <= x <= 1, this is the random variable for PCE estimation
% Transfer random x to the random values of uncertain parameters
[I_x,I_y,I_z,C_alpha,C_beta,SAP,R_0,k_z,delta_fL,delta_fR,V_wind,rho,Cd,FA,mu_rolling_1,mu_rolling_2,mu_dp,mu_lp,Vx,Vy,Vz,DOE_out] = rand_par_PCE(Q);

x1 = I_z;
x2 = C_alpha;
x3 = C_beta;
x4 = SAP;
x5 = mu_rolling_1;
x6 = mu_rolling_2;

% sub Collocation point in polynomial
H = matlabFunction(PHI');
parfor k = 1:length(x)
    phi(k,:) = H(x(k,1),x(k,2),x(k,3),x(k,4),x(k,5),x(k,6));
end
% Define polynomial coefficient of x by Least Square
INV = inv(phi'*phi)*phi';

C_x1 = double(INV*x1);
C_x2 = double(INV*x2);
C_x3 = double(INV*x3);
C_x4 = double(INV*x4);
C_x5 = double(INV*x5);
C_x6 = double(INV*x6);
toc

% Reconstruct model input in PCE form adding calculated coefficient
x1_re = PHI'*C_x1;
x2_re = PHI'*C_x2;
x3_re = PHI'*C_x3;
x4_re = PHI'*C_x4;
x5_re = PHI'*C_x5;
x6_re = PHI'*C_x6;

%% KALMAN FILTER
% ============================================================================
%% Real Model Parameters
terminal_time = 21;
desire_piece = 5;
k = terminal_time/desire_piece;
% Ream model Parameters
I_z_real = 48.8834;
C_alpha_real = 1413.4876;
C_beta_real = 56.6180;
SAP_real = 0.0029;
mu_rolling_1_real = 0.007;
mu_rolling_2_real = 5.2344e-04;

%% Closed-loop Drive
case_name = 'DLane';
sys_par = [7,0,0,1,1];
sim_time = terminal_time;
n_data_collect = 10000;
whole_par; % generate the nominal values of system parameters
% Drive once to generate maneuver command
[waypoint, require_velocity] = waypoints(sys_par,case_name); % create path following points
sim('tadpole_dynamic_close_loop'); % drive once, ouput its track's x and y points in cartisian coordinate
x_nominal = interp(linspace(0,length(x.Data),length(x.Data)),x.Data,n_data_collect);% x of each points
y_nominal = interp(linspace(0,length(y.Data),length(y.Data)),y.Data,n_data_collect);% y of each points
% Drive command
steer_command = steer_com.Data(1,:);
drive_command = drive_com.Data(1,:);
plot(x_nominal,y_nominal)

%% Opened-loop Drive
I_z = I_z_real;
C_alpha = C_alpha_real;
C_beta = C_beta_real;
SAP = SAP_real;
mu_rolling_1 = mu_rolling_1_real;
mu_rolling_2 = mu_rolling_2_real;
sim('tadpole_dynamic_open_loop2');
% Real Output Data
x_open_pt = interp(linspace(0,length(x_open.signals.values),length(x_open.signals.values)),x_open.signals.values,n_data_collect);
y_open_pt = interp(linspace(0,length(y_open.signals.values),length(y_open.signals.values)),y_open.signals.values,n_data_collect);
ax_open_pt = interp(linspace(0,length(a_x.signals.values),length(a_x.signals.values)),a_x.signals.values,n_data_collect);
ay_open_pt = interp(linspace(0,length(a_y.signals.values),length(a_y.signals.values)),a_y.signals.values,n_data_collect);
yaw_open_pt = interp(linspace(0,length(yaw_rate.signals.values),length(yaw_rate.signals.values)),yaw_rate.signals.values,n_data_collect);
time_pt = interp(linspace(0,length(x_open1.time),length(x_open1.time)),x_open1.time,n_data_collect);

H = eye(5);
R = zeros(5,5);
for j = 1:desire_piece
    tic
    % Simulation to desire time
    t_k = k*(j);
    % Create stochastic input; these are the parameters
    x1_kf = phi*C_x1;
    x2_kf = phi*C_x2;
    x3_kf = phi*C_x3;
    x4_kf = phi*C_x4;
    x5_kf = phi*C_x5;
    x6_kf = phi*C_x6;
    % Create corresponding simulation output
    sim_time = t_k;
    for p = 1:length(x1_kf)
        in(p) = Simulink.SimulationInput('tadpole_dynamic_open_loop2');%drive open loop vehicle
        in(p) = in(p).setVariable('l_1',l_1);
        in(p) = in(p).setVariable('l_2',l_2);
        in(p) = in(p).setVariable('w_1',w_1);
        in(p) = in(p).setVariable('w_2',w_2);
        in(p) = in(p).setVariable('h',h);
        in(p) = in(p).setVariable('p',p);
        in(p) = in(p).setVariable('m',m);
        in(p) = in(p).setVariable('g',9.81);
        in(p) = in(p).setVariable('tire_r',0.3302);
        in(p) = in(p).setVariable('tire_wdt',0.025);
        in(p) = in(p).setVariable('I_x',I_x);
        in(p) = in(p).setVariable('I_y',I_y);
        in(p) = in(p).setVariable('I_z',x1_kf(p)); %x1
        in(p) = in(p).setVariable('C_alpha',x2_kf(p)); %x2
        in(p) = in(p).setVariable('C_beta',x3_kf(p)); %x3
        in(p) = in(p).setVariable('SAP',x4_kf(p)); %x4
        in(p) = in(p).setVariable('R_0',R_0);
        in(p) = in(p).setVariable('k_z',k_z);
        in(p) = in(p).setVariable('delta_fL',delta_fL);
        in(p) = in(p).setVariable('delta_fR',delta_fR);
        in(p) = in(p).setVariable('V_wind',V_wind);
        in(p) = in(p).setVariable('rho',rho);
        in(p) = in(p).setVariable('Cd',Cd);
        in(p) = in(p).setVariable('FA',FA);
        in(p) = in(p).setVariable('mu_rolling_1',x5_kf(p)); %x5
        in(p) = in(p).setVariable('mu_rolling_2',x5_kf(p)); %x6
        in(p) = in(p).setVariable('mu_dp',mu_dp);
        in(p) = in(p).setVariable('mu_lp',mu_lp);
        in(p) = in(p).setVariable('Vx',0);
        in(p) = in(p).setVariable('Vy',0);
        in(p) = in(p).setVariable('Vz',0);
        in(p) = in(p).setVariable('I_f_wheel',0.3);
        in(p) = in(p).setVariable('I_r_wheel',0.5);
        in(p) = in(p).setVariable('drive_command',drive_command);
        in(p) = in(p).setVariable('steer_command',steer_command);
        in(p) = in(p).setVariable('sim_time',sim_time);
    end
    % set_param('tadpole_dynamic_open_loop1','AlgebraicLoopSolver','TrustRegion')
    out = parsim(in,'ShowSimulationManager','on','ShowProgress','on');%parellel computation
    
    y_ax = zeros(length(DOE_out),1);
    y_ay = zeros(length(DOE_out),1);
    y_x = zeros(length(DOE_out),1);
    y_y = zeros(length(DOE_out),1);
    y_yaw = zeros(length(DOE_out),1);
    % output vectors
    for i = 1:length(DOE_out)
        y_ax(i) = out(i).a_x.Data(end);
        y_ay(i) = out(i).a_y.Data(end);
        y_yaw(i) = out(i).yaw_rate.Data(end);
        y_x(i) = out(i).x_open.Data(end);
        y_y(i) = out(i).y_open.Data(end);
    end
    % real model output
    z_ax = interp1(a_x.Time,a_x.Data,t_k);
    z_ay = interp1(a_y.Time,a_y.Data,t_k);
    z_x = interp1(x_open.Time,x_open.Data,t_k);
    z_y = interp1(y_open.Time,y_open.Data,t_k);
    z_yaw = interp1(yaw_rate.Time,yaw_rate.Data,t_k);
    
    z_k = [z_ax;z_ay;z_yaw;z_x;z_y];
    
    R_ax = max(10^(-12),(0.03*z_ax)^2);
    R_ay = max(10^(-12),(0.03*z_ay)^2);
    R_yaw = max(10^(-12),((2.5/180*pi)^2));
    R_x = max(10^(-12),(0.025)^2);
    R_y = max(10^(-12),(0.025)^2);
    
    R = [R_ax 0 0 0 0;
         0   R_ay 0 0 0;
         0  0 R_yaw 0 0;
         0  0  0  R_x 0;
         0  0  0  0  R_y];
    
    % PCE of y_sim_kf
    C_y_ax_kf = double(inv(phi'*phi)*phi'*y_ax);
    C_y_ay_kf = double(inv(phi'*phi)*phi'*y_ay);
    C_y_x_kf = double(inv(phi'*phi)*phi'*y_x);
    C_y_y_kf = double(inv(phi'*phi)*phi'*y_y);
    C_y_yaw_kf = double(inv(phi'*phi)*phi'*y_yaw);
    
    
    
    % Calculate Covariance Matrix
    P_x1 = cov([y_ax,y_ay,y_yaw,y_x,y_y,x1_kf]); %6x6
    P_x2 = cov([y_ax,y_ay,y_yaw,y_x,y_y,x2_kf]); %6x6
    P_x3 = cov([y_ax,y_ay,y_yaw,y_x,y_y,x3_kf]); %6x6
    P_x4 = cov([y_ax,y_ay,y_yaw,y_x,y_y,x4_kf]); %6x6
    P_x5 = cov([y_ax,y_ay,y_yaw,y_x,y_y,x5_kf]); %6x6
    P_x6 = cov([y_ax,y_ay,y_yaw,y_x,y_y,x6_kf]); %6x6
    
    P_x1_y = P_x1(6,1:5); %1x5
    P_x2_y = P_x2(6,1:5); %1x5
    P_x3_y = P_x3(6,1:5); %1x5
    P_x4_y = P_x4(6,1:5); %1x5
    P_x5_y = P_x5(6,1:5); %1x5
    P_x6_y = P_x6(6,1:5); %1x5
    
    P_yy = P_x1(1:5,1:5);
    
    
    
    % Update coefficient of each polynomial at a time
    for i = 1:length(PHI)
        % Kernal Delta function
        if i ==1
            delta = 1;
        else
            delta = 0;
        end
        C_y_kf = [C_y_ax_kf(i);C_y_ay_kf(i);C_y_yaw_kf(i);C_y_x_kf(i);C_y_y_kf(i)];
        
        C_x1_a(i) = C_x1(i)+P_x1_y*H'*inv(R+H*P_yy*H)*(z_k*delta-H*C_y_kf);
        C_x2_a(i) = C_x2(i)+P_x2_y*H'*inv(R+H*P_yy*H)*(z_k*delta-H*C_y_kf);
        C_x3_a(i) = C_x3(i)+P_x3_y*H'*inv(R+H*P_yy*H)*(z_k*delta-H*C_y_kf);
        C_x4_a(i) = C_x4(i)+P_x4_y*H'*inv(R+H*P_yy*H)*(z_k*delta-H*C_y_kf);
        C_x5_a(i) = C_x5(i)+P_x5_y*H'*inv(R+H*P_yy*H)*(z_k*delta-H*C_y_kf);
        C_x6_a(i) = C_x6(i)+P_x6_y*H'*inv(R+H*P_yy*H)*(z_k*delta-H*C_y_kf);           
    end
        C_x1 = C_x1_a';
        C_x2 = C_x2_a';
        C_x3 = C_x3_a';
        C_x4 = C_x4_a';
        C_x5 = C_x5_a';
        C_x6 = C_x6_a';
    
       
end

% Final values
x1_final = phi*C_x1;
x2_final = phi*C_x2;
x3_final = phi*C_x3;
x4_final = phi*C_x4;
x5_final = phi*C_x5;
x6_final = phi*C_x6;




