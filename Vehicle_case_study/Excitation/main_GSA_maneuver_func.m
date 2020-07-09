%% SOBOL INDICES ON GLOBAL SENSITIVITY FOR TADPOLE THREE WHEEL VEHICLE
clc
close all

clear
tic
total_count = 0;
        %case 'inv_chirp'
            % A = [0.04,0.08,0.12,0.16,0.2]
            % w = [0.5,0.6,0.7,0.8,0.9,1]
            % gr = 0.1

for iii = 0.16:0.04:0.2
    for jjj = 0.5:0.1:1
        kkk = 0.01;
        total_count = total_count+1;
        %% 1. Complex Model Sampling
        % In this sampling stage, we first need to outcome the input and output
        % data as the learing data for kriging fitting.
        % input x: x is derived from sobol sequence depends on variable numbers and
        % sampling numbers.
        % output y: is the total squared error of output track and nominal values
        % set.
        
        runcycle = 1200; % total sampling times
        n_data_collect = 10000; % cut the output track into # pieces
        
        %---Nominal Set------------------------------------------------------
        % give maneuver parameters
        case_name = 'inv_chirp';
        sys_par = [iii,jjj,kkk]; %need a1, a2, b(better = 0), zoom_rate, zoom_length
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
        % insert parameters into simulink model
        % take the normalized number into design parameters of vehicle
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
        
        %% PCA
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
        
        
        % x
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
        
        % y
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
        
       
        GSA_N = [50000];
        for NUM = 1:length(GSA_N)
          
            %% Metamodel-Based Sensitivity Analysis
            N = GSA_N(NUM); % GSA monte carlo sampling numbers
            Q = sobolset(12,'skip',20000000);% create a different sobol sequence
            Q_ss = net(Q,N);
            P1 = Q_ss(:,1);
            P2 = Q_ss(:,2);
            P3 = Q_ss(:,3);
            P4 = Q_ss(:,4);
            P5 = Q_ss(:,5);
            P6 = Q_ss(:,6);
            Q1 = Q_ss(:,7);
            Q2 = Q_ss(:,8);
            Q3 = Q_ss(:,9);
            Q4 = Q_ss(:,10);
            Q5 = Q_ss(:,11);
            Q6 = Q_ss(:,12);
            
            P = [P1, P2, P3, P4, P5, P6];
            Q = [Q1, Q2, Q3, Q4, Q5, Q6];
            g_P = [];
            g_Q = [];
            
            % Sensitivity Analysis of Yaw Rate
            for i = 1:size(C_yaw,2)
                g_P = [];
                g_Q = [];
                for k = 1:length(P)/1000
                    g_P_temp = f_predictkrige(P((k-1)*1000+1:k*1000,:), kparam_yaw(i).par)';
                    g_Q_temp= f_predictkrige(Q((k-1)*1000+1:k*1000,:), kparam_yaw(i).par)';
                    g_P = [g_P g_P_temp];
                    g_Q = [g_Q g_Q_temp];
                end
                g_P_yaw(i).est = g_P;
                g_Q_yaw(i).est = g_Q;
                
                var_g = 1/N*(sum(g_P_yaw(i).est.^2))-(1/N*sum(g_P_yaw(i).est))^2; %Total variance
                
                for j = 1:6
                    PP = P;
                    QQ = Q;
                    PP(:,j) = QQ(:,j);
                    R(j) = struct('number',j,'in',PP); % will generate R matrix as sample matrix
                    g_RR = [];
                    for k = 1:length(PP)/1000
                        g_RR_temp = f_predictkrige(PP((k-1)*1000+1:k*1000,:), kparam_yaw(i).par)';
                        g_RR = [g_RR g_RR_temp];
                    end
                    g_R(j) = struct('no',j,'value',g_RR);% calulate each g(R^j)
                    C(j) = struct('no',j,'value',1/N*(sum(g_Q_yaw(i).est.*(g_R(j).value-g_P_yaw(i).est))));
                    D(j) = struct('no',j,'value',1/(2*N)*sum((g_P_yaw(i).est-g_R(j).value).^2));
                    Sj(j) = struct('main_effect',abs(C(j).value/var_g),'total_effect',abs(D(j).value/var_g));
                end
                S_index_yaw(i).index = Sj;
            end
            
            S_fusion_yaw = zeros(6,1);
            S_fusion_yaw_t = zeros(6,1);
            for i = 1:size(C_yaw,2)
                for j = 1:6
                    A(j,:) = (S_index_yaw(i).index(j).main_effect)*lambda_yaw(i);
                    B(j,:) = (S_index_yaw(i).index(j).total_effect)*lambda_yaw(i);
                end
                S_fusion_yaw = S_fusion_yaw + A/sum(lambda_yaw);
                S_fusion_yaw_t = S_fusion_yaw_t + B/sum(lambda_yaw);% Sensitivity of Data fusion
            end
            S_yaw_out = [S_fusion_yaw S_fusion_yaw_t];
            
            % Sensitivity Analysis of x
            for i = 1:size(C_x,2)
                g_P = [];
                g_Q = [];
                for k = 1:length(P)/1000
                    g_P_temp = f_predictkrige(P((k-1)*1000+1:k*1000,:), kparam_x(i).par)';
                    g_Q_temp= f_predictkrige(Q((k-1)*1000+1:k*1000,:), kparam_x(i).par)';
                    g_P = [g_P g_P_temp];
                    g_Q = [g_Q g_Q_temp];
                end
                g_P_x(i).est = g_P;
                g_Q_x(i).est = g_Q;
                
                var_g = 1/N*(sum(g_P_x(i).est.^2))-(1/N*sum(g_P_x(i).est))^2; %Total variance
                
                for j = 1:6
                    PP = P;
                    QQ = Q;
                    PP(:,j) = QQ(:,j);
                    R(j) = struct('number',j,'in',PP); % will generate R matrix as sample matrix
                    g_RR = [];
                    for k = 1:length(PP)/1000
                        g_RR_temp = f_predictkrige(PP((k-1)*1000+1:k*1000,:), kparam_x(i).par)';
                        g_RR = [g_RR g_RR_temp];
                    end
                    g_R(j) = struct('no',j,'value',g_RR);% calulate each g(R^j)
                    C(j) = struct('no',j,'value',1/N*(sum(g_Q_x(i).est.*(g_R(j).value-g_P_x(i).est))));
                    D(j) = struct('no',j,'value',1/(2*N)*sum((g_P_x(i).est-g_R(j).value).^2));
                    Sj(j) = struct('main_effect',abs(C(j).value/var_g),'total_effect',abs(D(j).value/var_g));
                end
                S_index_x(i).index = Sj;
            end
            
            S_fusion_x = zeros(6,1);
            S_fusion_x_t = zeros(6,1);
            for i = 1:size(C_x,2)
                for j = 1:6
                    A_x(j,:) = (S_index_x(i).index(j).main_effect)*lambda_x(i);
                    B_x(j,:) = (S_index_x(i).index(j).total_effect)*lambda_x(i);
                end
                S_fusion_x = S_fusion_x + A_x/sum(lambda_x); % Sensitivity of Data fusion
                S_fusion_x_t = S_fusion_x_t + B_x/sum(lambda_x);% Sensitivity of Data fusion
            end
            S_x_out = [S_fusion_x S_fusion_x_t];
            % Sensitivity Analysis of y
            for i = 1:size(C_y,2)
                g_P = [];
                g_Q = [];
                for k = 1:length(P)/1000
                    g_P_temp = f_predictkrige(P((k-1)*1000+1:k*1000,:), kparam_y(i).par)';
                    g_Q_temp= f_predictkrige(Q((k-1)*1000+1:k*1000,:), kparam_y(i).par)';
                    g_P = [g_P g_P_temp];
                    g_Q = [g_Q g_Q_temp];
                end
                g_P_y(i).est = g_P;
                g_Q_y(i).est = g_Q;
                
                var_g = 1/N*(sum(g_P_y(i).est.^2))-(1/N*sum(g_P_y(i).est))^2; %Total variance
                
                for j = 1:6
                    PP = P;
                    QQ = Q;
                    PP(:,j) = QQ(:,j);
                    R(j) = struct('number',j,'in',PP); % will generate R matrix as sample matrix
                    g_RR = [];
                    for k = 1:length(PP)/1000
                        g_RR_temp = f_predictkrige(PP((k-1)*1000+1:k*1000,:), kparam_y(i).par)';
                        g_RR = [g_RR g_RR_temp];
                    end
                    g_R(j) = struct('no',j,'value',g_RR);% calulate each g(R^j)
                    C(j) = struct('no',j,'value',1/N*(sum(g_Q_y(i).est.*(g_R(j).value-g_P_y(i).est))));
                    D(j) = struct('no',j,'value',1/(2*N)*sum((g_P_y(i).est-g_R(j).value).^2));
                    Sj(j) = struct('main_effect',abs(C(j).value/var_g),'total_effect',abs(D(j).value/var_g));
                end
                S_index_y(i).index = Sj;
            end
            
            S_fusion_y = zeros(6,1);
            S_fusion_y_t = zeros(6,1);
            for i = 1:size(C_y,2)
                for j = 1:6
                    A_y(j,:) = (S_index_y(i).index(j).main_effect)*lambda_y(i);
                    B_y(j,:) = (S_index_y(i).index(j).total_effect)*lambda_y(i);
                end
                S_fusion_y = S_fusion_y + A_y/sum(lambda_y); % Sensitivity of Data fusion
                S_fusion_y_t = S_fusion_y_t + B_y/sum(lambda_y);% Sensitivity of Data fusion
            end
            S_y_out = [S_fusion_y S_fusion_y_t];
            
            S_data.GSA_num = N;
            S_data.out_x = S_x_out;
            S_data.out_y = S_y_out;
            S_data.out_yaw = S_yaw_out;
            S_data.index_x = S_index_x;
            S_data.index_y = S_index_y;
            S_data.index_yaw = S_index_yaw;
            S_data.lambda_x = lambda_x;
            S_data.lambda_y = lambda_y;
            S_data.lambda_yaw = lambda_yaw;
            
            
        end
        toc
        Sim_data(total_count).par = [iii jjj];
        Sim_data(total_count).data = S_data;
    end
end
