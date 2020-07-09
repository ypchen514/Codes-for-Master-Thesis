clear
clc
close all


%% Math Model Validation
sys_par = [500,0.001,1e4;
    500,0.0001,1e4;
    500,0.00005,1e4;
    500,0.001,1e5;
    500,0.0001,1e5;
    500,0.00005,1e5;
    1000,0.001,1e4;
    1000,0.0001,1e4;
    1000,0.00005,1e4;
    1000,0.001,1e5;
    1000,0.0001,1e5;
    1000,0.00005,1e5;
    2000,0.001,1e4;
    2000,0.0001,1e4;
    2000,0.00005,1e4;
    2000,0.001,1e5;
    2000,0.0001,1e5;
    2000,0.00005,1e5];

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
    % Q_ss_1 = lhsdesign(N_s,4,'smooth','off');
    Q = sobolset(4,'skip',20000000);% create a different sobol sequence
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
    total_krig_sample(s) = length(output);
    EGO_time(s) = toc
    
    %% Cross validation Error Calculation
    Sample_Q = sobolset(4,'skip',100000);% create a different sobol sequence
    Sobol_sample = net(Sample_Q,N);
    % Calculate Accurate Value of g
    yval = g(alpha,p1,p2,p3,p4,Sobol_sample);
    mu_val(s) = mean(yval);
    var_val(s) = var(yval);
    % Calculate Value of same position from Kriging
    ykrig = f_predictkrige(Sobol_sample, kparam);
    mu_krig(s) = mean(ykrig);
    var_krig(s) = var(ykrig);
    % Error between math model and kriging
    err = yval'-ykrig;
    mu_err(s) = mean(err);
    var_err(s) = var(err);
    % Error Percentage
    per_err = err./yval';
    mu_per(s) = mean(per_err);
    var_per(s) = var(per_err);
    
    Data(s) = struct('Ini_sampling',sys_par(s,1),'Max_Var',sys_par(s,2),'Val_sampling',sys_par(s,3),'EGO_Num',total_krig_sample(s),...
        'math_y',yval,'mean_math',mu_val(s),'var_math',var_val(s),'krig_y',ykrig,'mean_krig',mu_krig(s),'var_krig',var_krig(s),...
        'Err',err,'mean_err',mu_err(s),'var_err',var_err(s));
    
    figure
    subplot(2,2,1)
    hist(yval,100);
    title(['$Y$ from Model with',num2str(sys_par(s,3)),'samples'],'Interpreter','latex');
    xlabel('Error Value','Interpreter','latex');
    ylabel('Cumulated number','Interpreter','latex');
    subplot(2,2,2)
    hist(ykrig,100);
    title(['$Y$ from Kriging with',num2str(sys_par(s,3)),'samples'],'Interpreter','latex');
    xlabel('Error Value','Interpreter','latex');
    ylabel('Cumulated number','Interpreter','latex');
    subplot(2,2,3)
    hist(err,100);
    title(['Estimation Error - $\sigma^2 \leq$',num2str(sys_par(s,2))],'Interpreter','latex');
    xlabel(['$\mu=$ ',num2str(mu_err(s)),', $\delta^2 =$ ',num2str(var_err(s))],'Interpreter','latex');
    ylabel('Cumulated number','Interpreter','latex');
    subplot(2,2,4)
    hist(per_err,100);
    title(['Error Percentage- $\sigma^2 \leq$',num2str(sys_par(s,2))],'Interpreter','latex');
    xlabel(['$\mu=$ ',num2str(mu_per(s)),', $\delta^2 =$ ',num2str(var_per(s))],'Interpreter','latex');
    ylabel('Cumulated number','Interpreter','latex');
end
