%% === Validation of Kriging Estimation ===
clear
clc
close all
tic
% runcycle_array = [400,600,800,1000,1200,1400,1600,1800,2000];
runcycle_array = [1800];
for r = 1:length(runcycle_array)
    %% Training
    sys_par_data = [0.1,0.75,0.1,1,1];
    sys_par = sys_par_data;
    iii = sys_par(1);
    jjj = sys_par(2);
    runcycle = runcycle_array(r); % total sampling times
    n_data_collect = 10000; % cut the output track into # pieces
    tic
    %---Nominal Set------------------------------------------------------
    % give maneuver parameters
    case_name = 'inv_chirp';
    sys_par = [iii,jjj,0.1,1,1]; %need a1, a2, b(better = 0), zoom_rate, zoom_length
    %         sys_par = [iii,jjj];
    sim_time = 30;
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
    fraction = 200;% each fraction has 200 samples
    x_open = [];
    y_open = [];
    ax_open = [];
    ay_open = [];
    yaw_open = [];
    for j = 1:(runcycle/fraction)
        for i = 1:fraction
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
            x_open_pt(i,:) = interp(linspace(0,length(out(1,i).x_open.Data),length(out(1,i).x_open.Data)),out(1,i).x_open.Data,n_data_collect);
            y_open_pt(i,:) = interp(linspace(0,length(out(1,i).y_open.Data),length(out(1,i).y_open.Data)),out(1,i).y_open.Data,n_data_collect);
            ax_open_pt(i,:) = interp(linspace(0,length(out(1,i).a_x.Data),length(out(1,i).a_x.Data)),out(1,i).a_x.Data,n_data_collect);
            ay_open_pt(i,:) = interp(linspace(0,length(out(1,i).a_y.Data),length(out(1,i).a_y.Data)),out(1,i).a_y.Data,n_data_collect);
            yaw_open_pt(i,:) = interp(linspace(0,length(out(1,i).yaw_rate.Data),length(out(1,i).yaw_rate.Data)),out(1,i).yaw_rate.Data,n_data_collect);
            %         plot(x_open_pt(i,:),y_open_pt(i,:));
            hold on
        end
        x_open = [x_open;x_open_pt];
        y_open = [y_open;y_open_pt];
        ax_open = [ax_open;ax_open_pt];
        ay_open = [ay_open;ay_open_pt];
        yaw_open = [yaw_open;yaw_open_pt];
    end
    
    % PCA
    
    % Yaw rate
    [C_yaw,trace_yaw,lambda_yaw,L_yaw,index_eig_val_yaw] = PCA(yaw_open);
    for j = 1:size(C_yaw,1)
        for i = 1:size(C_yaw,2)
            Y(i,:)=sqrt(lambda_yaw(i))*L_yaw(:,index_eig_val_yaw(i))*C_yaw(j,i);
        end
        yaw_sum(j,:) = sum(Y);
    end
    yaw_re = yaw_sum + mean(yaw_open);
    % fit kriging
    parfor i = 1:size(C_yaw,2)
        input = DOE; % 0.5~1.5
        output = C_yaw(:,i);
        lb = 0.5*ones(1,size(input,2));
        ub = 1.5*ones(1,size(input,2));
        fittype = 1;
        SCFtype = 1;
        max_variance = 1;
        kparam_yaw(i).par = f_variogram_fit(input, output, lb, ub);
    end
    
    % % Ax
    % [C_ax,trace_ax,lambda_ax,L_ax,index_eig_val_ax] = PCA(ax_open);
    % for j = 1:size(C_ax,1)
    %     for i = 1:size(C_ax,2)
    %         Y(i,:)=sqrt(lambda_ax(i))*L_ax(:,index_eig_val_ax(i))*C_ax(j,i);
    %     end
    %     ax_sum(j,:) = sum(Y);
    % end
    % ax_re = ax_sum + mean(ax_open);
    %
    % % fit kriging
    % parfor i = 1:size(C_ax,2)
    %     input = DOE;
    %     output = C_ax(:,i);
    %     lb = 0.5*ones(1,size(input,2));
    %     ub = 1.5*ones(1,size(input,2));
    %     fittype = 1;
    %     SCFtype = 1;
    %     max_variance = 1;
    %     kparam_ax(i).par = f_variogram_fit(input, output, lb, ub);
    % end
    %
    % % Ay
    % [C_ay,trace_ay,lambda_ay,L_ay,index_eig_val_ay] = PCA(ay_open);
    % for j = 1:size(C_ay,1)
    %     for i = 1:size(C_ay,2)
    %         Y(i,:)=sqrt(lambda_ay(i))*L_ay(:,index_eig_val_ay(i))*C_ay(j,i);
    %     end
    %     ay_sum(j,:) = sum(Y);
    % end
    % ay_re = ay_sum + mean(ay_open);
    %
    % % fit kriging
    % parfor i = 1:size(C_ay,2)
    %     input = DOE;
    %     output = C_ay(:,i);
    %     lb = 0.5*ones(1,size(input,2));
    %     ub = 1.5*ones(1,size(input,2));
    %     fittype = 1;
    %     SCFtype = 1;
    %     max_variance = 1;
    %     kparam_ay(i).par = f_variogram_fit(input, output, lb, ub);
    % end
    
    % X
    [C_x,trace_x,lambda_x,L_x,index_eig_val_x] = PCA(x_open);
    for j = 1:size(C_x,1)
        for i = 1:size(C_x,2)
            Y(i,:)=sqrt(lambda_x(i))*L_x(:,index_eig_val_x(i))*C_x(j,i);
        end
        x_sum(j,:) = sum(Y);
    end
    x_re = x_sum + mean(x_open);
    
    % fit kriging
    parfor i = 1:size(C_x,2)
        input = DOE;
        output = C_x(:,i);
        lb = 0.5*ones(1,size(input,2));
        ub = 1.5*ones(1,size(input,2));
        fittype = 1;
        SCFtype = 1;
        max_variance = 1;
        kparam_x(i).par = f_variogram_fit(input, output, lb, ub);
    end
    
    % Y
    [C_y,trace_y,lambda_y,L_y,index_eig_val_y] = PCA(y_open);
    for j = 1:size(C_y,1)
        for i = 1:size(C_y,2)
            Y(i,:)=sqrt(lambda_y(i))*L_y(:,index_eig_val_y(i))*C_y(j,i);
        end
        y_sum(j,:) = sum(Y);
    end
    y_re = y_sum + mean(y_open);
    
    % fit kriging
    parfor i = 1:size(C_y,2)
        input = DOE;
        output = C_y(:,i);
        lb = 0.5*ones(1,size(input,2));
        ub = 1.5*ones(1,size(input,2));
        fittype = 1;
        SCFtype = 1;
        max_variance = 1;
        kparam_y(i).par = f_variogram_fit(input, output, lb, ub);
    end
    
    
    
    
    %% Validation
    % Data from real model
    sys_par_data = [0.1,0.75,0.1,1,1];
    sys_par = sys_par_data;
    iii = sys_par(1);
    jjj = sys_par(2);
    runcycle = 1200; % total sampling times
    n_data_collect = 10000; % cut the output track into # pieces
    tic
    %---Nominal Set------------------------------------------------------
    % give maneuver parameters
    case_name = 'inv_chirp';
    sys_par = [iii,jjj,0.1,1,1]; %need a1, a2, b(better = 0), zoom_rate, zoom_length
    %         sys_par = [iii,jjj];
    sim_time = 30;
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
    [I_x,I_y,I_z,C_alpha,C_beta,SAP,R_0,k_z,delta_fL,delta_fR,V_wind,rho,Cd,FA,mu_rolling_1,mu_rolling_2,mu_dp,mu_lp,Vx,Vy,Vz,DOE] = rand_par2(runcycle);
    fraction = 200;% each fraction has 200 samples
    x_open = [];
    y_open = [];
    ax_open = [];
    ay_open = [];
    yaw_open = [];
    for j = 1:(runcycle/fraction)
        for i = 1:fraction
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
            x_open_pt(i,:) = interp(linspace(0,length(out(1,i).x_open.Data),length(out(1,i).x_open.Data)),out(1,i).x_open.Data,n_data_collect);
            y_open_pt(i,:) = interp(linspace(0,length(out(1,i).y_open.Data),length(out(1,i).y_open.Data)),out(1,i).y_open.Data,n_data_collect);
            ax_open_pt(i,:) = interp(linspace(0,length(out(1,i).a_x.Data),length(out(1,i).a_x.Data)),out(1,i).a_x.Data,n_data_collect);
            ay_open_pt(i,:) = interp(linspace(0,length(out(1,i).a_y.Data),length(out(1,i).a_y.Data)),out(1,i).a_y.Data,n_data_collect);
            yaw_open_pt(i,:) = interp(linspace(0,length(out(1,i).yaw_rate.Data),length(out(1,i).yaw_rate.Data)),out(1,i).yaw_rate.Data,n_data_collect);
            plot(x_open_pt(i,:),y_open_pt(i,:));
            hold on
        end
        x_open = [x_open;x_open_pt];
        y_open = [y_open;y_open_pt];
        ax_open = [ax_open;ax_open_pt];
        ay_open = [ay_open;ay_open_pt];
        yaw_open = [yaw_open;yaw_open_pt];
    end
    
    
    %% Real Coefficient from  Simulation
    [C_yaw,trace_yaw,lambda_yaw,L_yaw,index_eig_val_yaw] = PCA(yaw_open);
    % [C_ax,trace_ax,lambda_ax,L_ax,index_eig_val_ax] = PCA(ax_open);
    % [C_ay,trace_ay,lambda_ay,L_ay,index_eig_val_ay] = PCA(ay_open);
    [C_x,trace_x,lambda_x,L_x,index_eig_val_x] = PCA(x_open);
    [C_y,trace_y,lambda_y,L_y,index_eig_val_y] = PCA(y_open);
    
    % Reformation from PCA
    % yaw
    for j = 1:size(C_yaw,1)
        for i = 1:size(C_yaw,2)
            Y(i,:)=sqrt(lambda_yaw(i))*L_yaw(:,index_eig_val_yaw(i))*C_yaw(j,i);
        end
        yaw_sum_re(j,:) = sum(Y);
    end
    yaw_re = yaw_sum_re + mean(yaw_open);
    
    % x
    for j = 1:size(C_x,1)
        for i = 1:size(C_x,2)
            Y(i,:)=sqrt(lambda_x(i))*L_x(:,index_eig_val_x(i))*C_x(j,i);
        end
        x_sum_re(j,:) = sum(Y);
    end
    x_re = x_sum_re + mean(x_open);
    
    % y
    for j = 1:size(C_y,1)
        for i = 1:size(C_y,2)
            Y(i,:)=sqrt(lambda_y(i))*L_y(:,index_eig_val_y(i))*C_y(j,i);
        end
        y_sum_re(j,:) = sum(Y);
    end
    y_re = y_sum_re + mean(y_open);
    
    
    
    
    % Data Estimate from Kriging
    for i = 1:length(lambda_yaw)
        C_yaw_krig(i).par = f_predictkrige(DOE, kparam_yaw(i).par);
        R_sq(i) = 1-(sum((C_yaw_krig(i).par-C_yaw(:,i)).^2))/(sum((C_yaw(:,i)-mean(C_yaw(:,i))).^2));
        if R_sq(i) < -1
            C_yaw_krig(i).par = -C_yaw_krig(i).par;
        end
        error_yaw(i).per = abs((C_yaw(:,i)-C_yaw_krig(i).par)./C_yaw(:,i))*100;
        R_sq_yaw(i) = 1-(sum((C_yaw_krig(i).par-C_yaw(:,i)).^2))/(sum((C_yaw(:,i)-mean(C_yaw(:,i))).^2));
        RAAE_yaw(i) = (sum(abs(C_yaw(:,i)-C_yaw_krig(i).par)))/(sqrt(length(C_yaw)*sum((C_yaw(:,i)-mean(C_yaw(:,i))).^2)));
    end
    
    R_sq_yaw_fusion = R_sq_yaw*lambda_yaw/sum(lambda_yaw);
    RAAE_yaw_fusion = RAAE_yaw*lambda_yaw/sum(lambda_yaw);
    
    % for i = 1:length(lambda_ax)
    % C_ax_krig(i).par = f_predictkrige(DOE, kparam_ax(i).par);
    %  R_sq(i) = 1-(sum((C_ax_krig(i).par-C_ax(:,i)).^2))/(sum((C_ax(:,i)-mean(C_ax(:,i))).^2));
    %     if R_sq(i) < 0
    %         C_ax_krig(i).par = -C_ax_krig(i).par;
    %     end
    % error_ax(i).per = abs((C_ax(:,i)-C_ax_krig(i).par)./C_ax(:,i))*100;
    % R_sq_ax(i) = 1-(sum((C_ax_krig(i).par-C_ax(:,i)).^2))/(sum((C_ax(:,i)-mean(C_ax(:,i))).^2));
    % RAAE_ax(i) = (sum(abs(C_ax(:,i)-C_ax_krig(i).par)))/(sqrt(length(C_ax)*sum((C_ax(:,i)-mean(C_ax(:,i))).^2)));
    % end
    %
    % for i = 1:length(lambda_ay)
    % C_ay_krig(i).par = f_predictkrige(DOE, kparam_ay(i).par);
    % R_sq(i) = 1-(sum((C_ay_krig(i).par-C_ay(:,i)).^2))/(sum((C_ay(:,i)-mean(C_ay(:,i))).^2));
    %     if R_sq(i) < 0
    %         C_ay_krig(i).par = -C_ay_krig(i).par;
    %     end
    % error_ay(i).per = abs((C_ay(:,i)-C_ay_krig(i).par)./C_ay(:,i))*100;
    % R_sq_ay(i) = 1-(sum((C_ay_krig(i).par-C_ay(:,i)).^2))/(sum((C_ay(:,i)-mean(C_ay(:,i))).^2));
    % RAAE_ay(i) = (sum(abs(C_ay(:,i)-C_ay_krig(i).par)))/(sqrt(length(C_ay)*sum((C_ay(:,i)-mean(C_ay(:,i))).^2)));
    % end
    
    % X
    for i = 1:length(lambda_x)
        C_x_krig(i).par = f_predictkrige(DOE, kparam_x(i).par);
        R_sq(i) = 1-(sum((C_x_krig(i).par-C_x(:,i)).^2))/(sum((C_x(:,i)-mean(C_x(:,i))).^2));
        if R_sq(i) < 0
            C_x_krig(i).par = -C_x_krig(i).par;
        end
        error_x(i).per = abs((C_x(:,i)-C_x_krig(i).par)./(C_x(:,i)))*100;
        R_sq_x(i) = 1-(sum((C_x_krig(i).par-C_x(:,i)).^2))/(sum((C_x(:,i)-mean(C_x(:,i))).^2));
        RAAE_x(i) = (sum(abs(C_x(:,i)-C_x_krig(i).par)))/(sqrt(length(C_x)*sum((C_x(:,i)-mean(C_x(:,i))).^2)));
    end
    
    R_sq_x_fusion = R_sq_x*lambda_x/sum(lambda_x);
    RAAE_x_fusion = RAAE_x*lambda_x/sum(lambda_x);
    
    % Y
    for i = 1:length(lambda_y)
        C_y_krig(i).par = f_predictkrige(DOE, kparam_y(i).par);
        R_sq(i) = 1-(sum((C_y_krig(i).par-C_y(:,i)).^2))/(sum((C_y(:,i)-mean(C_y(:,i))).^2));
        if R_sq(i) < 0
            C_y_krig(i).par = -C_y_krig(i).par;
        end
        error_y(i).per = abs((C_y(:,i)-C_y_krig(i).par)./C_y(:,i))*100;
        R_sq_y(i) = 1-(sum((C_y_krig(i).par-C_y(:,i)).^2))/(sum((C_y(:,i)-mean(C_y(:,i))).^2));
        RAAE_y(i) = (sum(abs(C_y(:,i)-C_y_krig(i).par)))/(sqrt(length(C_y)*sum((C_y(:,i)-mean(C_y(:,i))).^2)));
    end
    
    R_sq_y_fusion = R_sq_y*lambda_y/sum(lambda_y);
    RAAE_y_fusion = RAAE_y*lambda_y/sum(lambda_y);
    
    Yaw_acc(r).test_N = runcycle;
    Yaw_acc(r).err_dist = error_yaw;
    Yaw_acc(r).RAAE = RAAE_yaw;
    Yaw_acc(r).RAAE_fusion = RAAE_yaw_fusion;
    Yaw_acc(r).R_sq = R_sq_yaw;
    Yaw_acc(r).R_sq_fusion = R_sq_yaw_fusion;
    
    x_acc(r).test_N = runcycle;
    x_acc(r).err_dist = error_x;
    x_acc(r).RAAE = RAAE_x;
    x_acc(r).RAAE_fusion = RAAE_x_fusion;
    x_acc(r).R_sq = R_sq_x;
    x_acc(r).R_sq_fusion = R_sq_x_fusion;
    
    y_acc(r).test_N = runcycle;
    y_acc(r).err_dist = error_y;
    y_acc(r).RAAE = RAAE_y;
    y_acc(r).RAAE_fusion = RAAE_y_fusion;
    y_acc(r).R_sq = R_sq_y;
    y_acc(r).R_sq_fusion = R_sq_y_fusion;
    
    
    
end

%% Reformation from Kriging
% yaw
for j = 1:size(C_yaw_krig(1).par,1)
    for i = 1:size(C_yaw_krig,2)
        Y(i,:)=sqrt(lambda_yaw(i))*L_yaw(:,index_eig_val_yaw(i))*C_yaw_krig(i).par(j);
    end
    yaw_sum_re_krig(j,:) = sum(Y);
end
yaw_re_krig = yaw_sum_re_krig + mean(yaw_open);

% x
for j = 1:size(C_x_krig(1).par,1)
    for i = 1:size(C_x_krig,2)
        Y(i,:)=sqrt(lambda_x(i))*L_x(:,index_eig_val_x(i))*C_x_krig(i).par(j);
    end
    x_sum_re_krig(j,:) = sum(Y);
end
x_re_krig = x_sum_re_krig + mean(x_open);

% y
for j = 1:size(C_y_krig(1).par,1)
    for i = 1:size(C_y_krig,2)
        Y(i,:)=sqrt(lambda_y(i))*L_y(:,index_eig_val_y(i))*C_y_krig(i).par(j);
    end
    y_sum_re_krig(j,:) = sum(Y);
end
y_re_krig = y_sum_re_krig + mean(y_open);

% Randomly select one set and plot
N = randi([0,runcycle]);
figure
subplot(1,3,1)
plot(x_open(N,:));
hold on
plot(x_re(N,:));
plot(x_re_krig(N,:));
xlabel('Data point','interpreter','latex')
ylabel('x location (meter)','interpreter','latex');
title('Comparison of $X$','Interpreter','latex');
legend('Real Data','Reformation from PCA','Reformation from Kriging','Interpreter','latex')

hold off

subplot(1,3,2)
plot(y_open(N,:));
hold on
plot(y_re(N,:));
plot(y_re_krig(N,:));
xlabel('Data point','interpreter','latex')
ylabel('y location (meter)','interpreter','latex');
title('Comparison of $Y$','Interpreter','latex');
legend('Real Data','Reformation from PCA','Reformation from Kriging','Interpreter','latex')
hold off


subplot(1,3,3)
plot(yaw_open(N,:));
hold on
plot(yaw_re(N,:));
plot(yaw_re_krig(N,:));
xlabel('Data point','interpreter','latex')
ylabel('yaw rate (rad/s)','interpreter','latex');
title('Comparison of yaw rate','Interpreter','latex');
legend('Real Data','Reformation from PCA','Reformation from Kriging','Interpreter','latex')
hold off




i = 1;
% hist x
figure
subplot(1,4,1)
A = sort(x_acc(i).err_dist(1).per);
hist(A(1:0.9*runcycle),20);
xlabel('Error Percentage','interpreter','latex')
ylabel('Cumulated Number','interpreter','latex')
title('Error distribution of Kriging on $C_1$','interpreter','latex');
subplot(1,4,2)
A = sort(x_acc(i).err_dist(2).per);
hist(A(1:0.9*runcycle),20);
xlabel('Error Percentage','interpreter','latex')
ylabel('Cumulated Number','interpreter','latex')
title('Error distribution of Kriging on $C_2$','interpreter','latex');
subplot(1,4,3)
A = sort(x_acc(i).err_dist(3).per);
hist(A(1:0.9*runcycle),20);
xlabel('Error Percentage','interpreter','latex')
ylabel('Cumulated Number','interpreter','latex')
title('Error distribution of Kriging on $C_3$','interpreter','latex');
subplot(1,4,4)
A = sort(x_acc(i).err_dist(4).per);
hist(A(1:0.9*runcycle),20);
xlabel('Error Percentage','interpreter','latex')
ylabel('Cumulated Number','interpreter','latex')
title('Error distribution of Kriging on $C_4$','interpreter','latex');

% hist y
figure
subplot(1,3,1)
A = sort(y_acc(i).err_dist(1).per);
hist(A(1:0.9*runcycle),20);
xlabel('Error Percentage','interpreter','latex')
ylabel('Cumulated Number','interpreter','latex')
title('Error distribution of Kriging on $C_1$','interpreter','latex');
subplot(1,3,2)
A = sort(y_acc(i).err_dist(2).per);
hist(A(1:0.9*runcycle),20);
xlabel('Error Percentage','interpreter','latex')
ylabel('Cumulated Number','interpreter','latex')
title('Error distribution of Kriging on $C_2$','interpreter','latex');
subplot(1,3,3)
A = sort(y_acc(i).err_dist(3).per);
hist(A(1:0.9*runcycle),20);
xlabel('Error Percentage','interpreter','latex')
ylabel('Cumulated Number','interpreter','latex')
title('Error distribution of Kriging on $C_3$','interpreter','latex');

% yaw estimate
figure
subplot(1,3,1)
A = sort(Yaw_acc(i).err_dist(1).per);
hist(A(1:0.9*runcycle),20);
xlabel('Error Percentage','interpreter','latex')
ylabel('Cumulated Number','interpreter','latex')
title('Error distribution of Kriging on $C_1$','interpreter','latex');
subplot(1,3,2)
A = sort(Yaw_acc(i).err_dist(2).per);
hist(A(1:0.9*runcycle),20);
xlabel('Error Percentage','interpreter','latex')
ylabel('Cumulated Number','interpreter','latex')
title('Error distribution of Kriging on $C_2$','interpreter','latex');
subplot(1,3,3)
A = sort(Yaw_acc(i).err_dist(3).per);
hist(A(1:0.9*runcycle),20);
xlabel('Error Percentage','interpreter','latex')
ylabel('Cumulated Number','interpreter','latex')
title('Error distribution of Kriging on $C_3$','interpreter','latex');
% subplot(1,4,4)
% A = sort(Yaw_acc(i).err_dist(4).per);
% hist(A(1:0.9*runcycle),20);
% xlabel('Error Percentage','interpreter','latex')
% ylabel('Cumulated Number','interpreter','latex')
% title('Error distribution of Kriging on $C_4$','interpreter','latex');

%% Reformation dataset comparison
% dataset of x
figure
subplot(1,3,3)
for i = 1:runcycle
    plot(x_re_krig(i,:));
    hold on
end
title('Reformed Data from Kriged-PCA','interpreter','latex','fontsize',14)
xlabel('Data point','interpreter','latex');
ylabel('$X$ location (m)','interpreter','latex');

subplot(1,3,2)
for i = 1:runcycle
    plot(x_re(i,:));
    hold on
end
title('Reformed Data from PCA','interpreter','latex','fontsize',14)
xlabel('Data point','interpreter','latex');
ylabel('$X$ location (m)','interpreter','latex');

subplot(1,3,1)
for i = 1:runcycle
    plot(x_open(i,:));
    hold on
end
title('Original Data','interpreter','latex','fontsize',14)
xlabel('Data point','interpreter','latex');
ylabel('$X$ location (m)','interpreter','latex');


% dataset of y
figure
subplot(1,3,3)
for i = 1:runcycle
    plot(y_re_krig(i,:));
    hold on
end
title('Reformed Data from Kriged-PCA','interpreter','latex','fontsize',14)
xlabel('Data point','interpreter','latex');
ylabel('$Y$ location (m)','interpreter','latex');

subplot(1,3,2)
for i = 1:runcycle
    plot(y_re(i,:));
    hold on
end
title('Reformed Data from PCA','interpreter','latex','fontsize',14)
xlabel('Data point','interpreter','latex');
ylabel('$Y$ location (m)','interpreter','latex');

subplot(1,3,1)
for i = 1:runcycle
    plot(y_open(i,:));
    hold on
end
title('Original Data','interpreter','latex','fontsize',14)
xlabel('Data point','interpreter','latex');
ylabel('$Y$ location (m)','interpreter','latex');



% dataset of yaw
figure
subplot(1,3,3)
for i = 1:runcycle
    plot(yaw_re_krig(i,:));
    hold on
end
title('Reformed Data from Kriged-PCA','interpreter','latex','fontsize',14)
xlabel('Data point','interpreter','latex');
ylabel('Yaw rate (rad/s)','interpreter','latex');

subplot(1,3,2)
for i = 1:runcycle
    plot(yaw_re(i,:));
    hold on
end
title('Reformed Data from PCA','interpreter','latex','fontsize',14)
xlabel('Data point','interpreter','latex');
ylabel('Yaw rate (rad/s)','interpreter','latex');

subplot(1,3,1)
for i = 1:runcycle
    plot(yaw_open(i,:));
    hold on
end
title('Original Data','interpreter','latex','fontsize',14)
xlabel('Data point','interpreter','latex');
ylabel('Yaw rate (rad/s)','interpreter','latex');


toc

% x C plot
figure
subplot(4,1,1)
hold on
plot(C_x(1:200,1),'o')
plot(C_x_krig(1).par(1:200),'r.','markersize',12)
plot(C_x(1:200,1))
hold off
legend('$C_1$ from Calculation','$C_1$ from Kriging Prediction','Interpreter','latex');
xlabel('Data points','Interpreter','latex');
ylabel('Coefficient value','interpreter','latex');
title('Fitting Accuracy of $C_1$','interpreter','latex','fontsize',14);
subplot(4,1,2)
hold on
plot(C_x(1:200,2),'o')
plot(C_x_krig(2).par(1:200),'r.','markersize',12)
plot(C_x(1:200,2))
hold off
legend('$C_2$ from Calculation','$C_2$ from Kriging Prediction','Interpreter','latex');
xlabel('Data points','Interpreter','latex');
ylabel('Coefficient value','interpreter','latex');
title('Fitting Accuracy of $C_2$','interpreter','latex','fontsize',14);
subplot(4,1,3)
hold on
plot(C_x(1:200,3),'o')
plot(C_x_krig(3).par(1:200),'r.','markersize',12)
plot(C_x(1:200,3))
hold off
legend('$C_3$ from Calculation','$C_3$ from Kriging Prediction','Interpreter','latex');
xlabel('Data points','Interpreter','latex');
ylabel('Coefficient value','interpreter','latex');
title('Fitting Accuracy of $C_3$','interpreter','latex','fontsize',14);
subplot(4,1,4)
hold on
plot(C_x(1:200,4),'o')
plot(C_x_krig(4).par(1:200),'r.','markersize',12)
plot(C_x(1:200,4))
hold off
legend('$C_4$ from Calculation','$C_4$ from Kriging Prediction','Interpreter','latex');
xlabel('Data points','Interpreter','latex');
ylabel('Coefficient value','interpreter','latex');
title('Fitting Accuracy of $C_4$','interpreter','latex','fontsize',14);


% y C plot
figure
subplot(3,1,1)
hold on
plot(C_y(1:200,1),'o')
plot(C_y_krig(1).par(1:200),'r.','markersize',12)
plot(C_y(1:200,1))
hold off
legend('$C_1$ from Calculation','$C_1$ from Kriging Prediction','Interpreter','latex');
xlabel('Data points','Interpreter','latex');
ylabel('Coefficient value','interpreter','latex');
title('Fitting Accuracy of $C_1$','interpreter','latex','fontsize',14);
subplot(3,1,2)
hold on
plot(C_y(1:200,2),'o')
plot(C_y_krig(2).par(1:200),'r.','markersize',12)
plot(C_y(1:200,2))
hold off
legend('$C_2$ from Calculation','$C_2$ from Kriging Prediction','Interpreter','latex');
xlabel('Data points','Interpreter','latex');
ylabel('Coefficient value','interpreter','latex');
title('Fitting Accuracy of $C_2$','interpreter','latex','fontsize',14);
subplot(3,1,3)
hold on
plot(C_y(1:200,3),'o')
plot(C_y_krig(3).par(1:200),'r.','markersize',12)
plot(C_y(1:200,3))
hold off
legend('$C_3$ from Calculation','$C_3$ from Kriging Prediction','Interpreter','latex');
xlabel('Data points','Interpreter','latex');
ylabel('Coefficient value','interpreter','latex');
title('Fitting Accuracy of $C_3$','interpreter','latex','fontsize',14);


% yaw C plot
figure
subplot(3,1,1)
hold on
plot(C_yaw(1:200,1),'o')
plot(C_yaw_krig(1).par(1:200),'r.','markersize',12)
plot(C_yaw(1:200,1))
hold off
legend('$C_1$ from Calculation','$C_1$ from Kriging Prediction','Interpreter','latex');
xlabel('Data points','Interpreter','latex');
ylabel('Coefficient value','interpreter','latex');
title('Fitting Accuracy of $C_1$','interpreter','latex','fontsize',14);
subplot(3,1,2)
hold on
plot(C_yaw(1:200,2),'o')
plot(C_yaw_krig(2).par(1:200),'r.','markersize',12)
plot(C_yaw(1:200,2))
hold off
legend('$C_2$ from Calculation','$C_2$ from Kriging Prediction','Interpreter','latex');
xlabel('Data points','Interpreter','latex');
ylabel('Coefficient value','interpreter','latex');
title('Fitting Accuracy of $C_2$','interpreter','latex','fontsize',14);
subplot(3,1,3)
hold on
plot(C_yaw(1:200,3),'o')
plot(C_yaw_krig(3).par(1:200),'r.','markersize',12)
plot(C_yaw(1:200,3))
hold off
legend('$C_3$ from Calculation','$C_3$ from Kriging Prediction','Interpreter','latex');
xlabel('Data points','Interpreter','latex');
ylabel('Coefficient value','interpreter','latex');
title('Fitting Accuracy of $C_3$','interpreter','latex','fontsize',14);
% subplot(4,1,4)
% hold on
% plot(C_yaw(1:200,4),'o')
% plot(C_yaw_krig(4).par(1:200),'r.','markersize',12)
% plot(C_yaw(1:200,4))
% hold off
% legend('$C_4$ from Calculation','$C_4$ from Kriging Prediction','Interpreter','latex');
% xlabel('Data points','Interpreter','latex');
% ylabel('Coefficient value','interpreter','latex');
% title('Fitting Accuracy of $C_4$','interpreter','latex','fontsize',14);
