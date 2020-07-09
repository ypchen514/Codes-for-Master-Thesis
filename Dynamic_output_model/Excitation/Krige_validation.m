%% Validation of Kriged-PCA Based Sensitivity Analysis
clc
clear
close all

tic
N_sample_e = 200;
%% Model parameters
% Operation parameters
x1 = 10;
x2 = 0.5;
x3 = 10;
% Random Parameters by Low-discrepancy Sampling
Q = sobolset(8,'skip',200000);% create a different sobol sequence
Q_ss = net(Q,N_sample_e);
P_a = (Q_ss(:,1)-0.5)*0.1+1; % uncertainty of a
P_b = (Q_ss(:,2)-0.5)*0.1+1; % uncertainty of b
P_c = (Q_ss(:,3)-0.5)*0.1+1; % uncertainty of c
P_n = norminv(Q_ss(:,4),0,0.01); % sensing error e_e


Q_a = (Q_ss(:,4)-0.5)*0.1+1; % uncertainty
Q_b = (Q_ss(:,5)-0.5)*0.1+1; % uncertainty
Q_c = (Q_ss(:,6)-0.5)*0.1+1; % uncertainty
Q_n = norminv(Q_ss(:,8),0,0.01); % sensing error e_e

%% Build Kriging Replacing the Math System
% Kriging input: P_a,P_b,P_c,P_n
% Kriging output: C1,C2,C3,C4,C5
Krig_in = [Q_ss(:,1) Q_ss(:,2) Q_ss(:,3) Q_ss(:,4)]; % Q_ss from 0~1
[C_1_out, C_2_out, C_3_out, C_4_out,C_5_out,lambda_g,y_e] = g(P_a,P_b,P_c,P_n,x1,x2,x3); %P varies as its possible range
plot(y_e(:,:)');

% Build Kriging
lb = [0,0,0,0];
ub = [1,1,1,1];
fittype = 1;
SCFtype = 1;
max_variance = 1;
run_num = 0;

kparam_C1 = f_variogram_fit(Krig_in, C_1_out, lb, ub);
kparam_C2 = f_variogram_fit(Krig_in, C_2_out, lb, ub);
kparam_C3 = f_variogram_fit(Krig_in, C_3_out, lb, ub);
kparam_C4 = f_variogram_fit(Krig_in, C_4_out, lb, ub);
kparam_C5 = f_variogram_fit(Krig_in, C_5_out, lb, ub);

%% EGO
CVar = 0.05; % converge variance  = 0.05
trial = 0;
[x_next, variance_next] = f_find_kriging_max_variance(kparam_C1);
while abs(variance_next) >= CVar && trial <= 50
    Krig_in = [Krig_in;[x_next(1) x_next(2) x_next(3) x_next(4)]]; 
    P_a = [P_a;(x_next(1)-0.5)*0.1+1];
    P_b = [P_b;(x_next(2)-0.5)*0.1+1];
    P_c = [P_c;(x_next(3)-0.5)*0.1+1];
    P_n = [P_n;normrnd(0,0.01)];
    [C_1_out, C_2_out, C_3_out, C_4_out,C_5_out,lambda_g,y_e] = g(P_a,P_b,P_c,P_n,x1,x2,x3); %P varies as its possible range
    % Build Kriging
    lb = [0,0,0,0];
    ub = [1,1,1,1];
    fittype = 1;
    SCFtype = 1;
    max_variance = 1;
    run_num = 0;
    
    kparam_C1 = f_variogram_fit(Krig_in, C_1_out, lb, ub);
    kparam_C2 = f_variogram_fit(Krig_in, C_2_out, lb, ub);
    kparam_C3 = f_variogram_fit(Krig_in, C_3_out, lb, ub);
    kparam_C4 = f_variogram_fit(Krig_in, C_4_out, lb, ub);
    kparam_C5 = f_variogram_fit(Krig_in, C_5_out, lb, ub);
    [x_next, variance_next] = f_find_kriging_max_variance(kparam_C1);
    trial = trial + 1;
end



toc1 = toc

N_krig_sample = 1000;
Q = sobolset(4,'skip',4000000);% create a different sobol sequence
Q_krig = net(Q,N_krig_sample);

P_a = (Q_krig(:,1)-0.5)*0.1+1; % uncertainty of a
P_b = (Q_krig(:,2)-0.5)*0.1+1; % uncertainty of b
P_c = (Q_krig(:,3)-0.5)*0.1+1; % uncertainty of c
P_n = norminv(Q_krig(:,4),0,0.01); % sensing error e_e

[C_1_out, C_2_out, C_3_out, C_4_out,C_5_out,lambda,y_e] = g(P_a,P_b,P_c,P_n,x1,x2,x3);% Coefficient calculated from real model  

Yc = y_e-mean(y_e);
T = cov(Yc); % Covariance Matrix
Tot_Iner = trace(T); % Total inertial = trace of covariance matrix
[L eig_val] = eig(T);
eig_val = eig(T);
[lambda index_eig_val] = sort(eig_val,'descend');

PC_p1 = f_predictkrige(Q_krig(:,1:4),kparam_C1);
PC_p2 = f_predictkrige(Q_krig(:,1:4),kparam_C2);
PC_p3 = f_predictkrige(Q_krig(:,1:4),kparam_C3);
PC_p4 = f_predictkrige(Q_krig(:,1:4),kparam_C4);
PC_p5 = f_predictkrige(Q_krig(:,1:4),kparam_C5);

%% Time-dependent output validation
% Output from real model
C = [C_1_out C_2_out C_3_out C_4_out C_5_out];
for i = 1:5
   Y(i,:)=sqrt(lambda(i))*L(:,index_eig_val(i))*C(1,i); 
end    
Y_re_exp = sum(Y)+mean(y_e);

% Output from Kriging
x_input = [P_a P_b P_c P_n];
PC = [PC_p1 PC_p2 PC_p3 PC_p4 PC_p5];
for i = 1:5
   Y_krig(i,:)=sqrt(lambda(i))*L(:,index_eig_val(i))*PC(1,i); 
end    
Y_re_krig = sum(Y_krig)+mean(y_e);
t = linspace(0,15,1500);

figure
subplot(1,2,1)
hold on
plot(t,y_e(1,:),'linewidth',3);
plot(t,Y_re_exp,'r','linewidth',1.5);
plot(t,Y_re_krig,'b--','linewidth',1.5);
hold off
xlabel('Time (s)', 'interpreter', 'latex');
ylabel('$y$','interpreter','latex');
legend('Simulation Output','Reconstruct from $C_i$','Reconstruct from $\hat{C}_i$','interpreter','latex')
title('Trajectory Reformation and Comparison','interpreter','latex');



subplot(1,2,2)
hold on
plot(t,y_e(1,:),'linewidth',3);
plot(t,Y_re_exp,'r','linewidth',1.5);
plot(t,Y_re_krig,'b--','linewidth',1.5);
hold off
xlabel('Time (s)', 'interpreter', 'latex');
ylabel('$y$','interpreter','latex');
legend('Simulation Output','Reconstruct from $C_i$','Reconstruct from $\hat{C}_i$','interpreter','latex')
title('Trajectory Reformation and Comparison (Partial)','interpreter','latex');

% R-square
R2_C1 = 1-(sum((PC_p1-C_1_out).^2))/(sum((C_1_out-mean(C_1_out)).^2));
R2_C2 = 1-(sum((PC_p2-C_2_out).^2))/(sum((C_2_out-mean(C_2_out)).^2));
R2_C3 = 1-(sum((PC_p3-C_3_out).^2))/(sum((C_3_out-mean(C_3_out)).^2));
R2_C4 = 1-(sum((PC_p4-C_4_out).^2))/(sum((C_4_out-mean(C_4_out)).^2));
R2_C5 = 1-(sum((PC_p5-C_5_out).^2))/(sum((C_5_out-mean(C_5_out)).^2));

% RAAE
RAAE_C1 = (sum(abs(C_1_out-PC_p1)))/(sqrt(N_krig_sample*sum((C_1_out-mean(C_1_out)).^2)));
RAAE_C2 = (sum(abs(C_2_out-PC_p2)))/(sqrt(N_krig_sample*sum((C_2_out-mean(C_2_out)).^2)));
RAAE_C3 = (sum(abs(C_3_out-PC_p3)))/(sqrt(N_krig_sample*sum((C_3_out-mean(C_3_out)).^2)));
RAAE_C4 = (sum(abs(C_4_out-PC_p4)))/(sqrt(N_krig_sample*sum((C_4_out-mean(C_4_out)).^2)));
RAAE_C5 = (sum(abs(C_5_out-PC_p5)))/(sqrt(N_krig_sample*sum((C_5_out-mean(C_5_out)).^2)));

% Error plot
err_C1 = abs((PC_p1 - C_1_out)./C_1_out*100);
err_C2 = abs((PC_p2 - C_2_out)./C_2_out*100);
err_C3 = abs((PC_p3 - C_3_out)./C_3_out*100);
err_C4 = abs((PC_p4 - C_4_out)./C_4_out*100);
err_C5 = abs((PC_p5 - C_5_out)./C_5_out*100);

figure
subplot(2,3,1)
hist(err_C1,2000);
title('Error Percentage of $C_1$','Interpreter','latex','Fontsize',18);
xlabel('Error Percentage ','Interpreter','latex','Fontsize',16);

subplot(2,3,2);
hist(err_C2,500);
title('Error Percentage of $C_2$','Interpreter','latex','Fontsize',18);
xlabel('Error Percentage ','Interpreter','latex','Fontsize',16);

subplot(2,3,3);
hist(err_C3,5000);
title('Error Percentage of $C_3$','Interpreter','latex','Fontsize',18);
xlabel('Error Percentage ','Interpreter','latex','Fontsize',16);

subplot(2,3,4);
hist(err_C4,10000);
title('Error Percentage of $C_4$','Interpreter','latex','Fontsize',18);
xlabel('Error Percentage ','Interpreter','latex','Fontsize',16);

subplot(2,3,5);
hist(err_C5,2000);
title('Error Percentage of $C_5$','Interpreter','latex','Fontsize',18);
xlabel('Error Percentage ','Interpreter','latex','Fontsize',16);

output = [N_sample_e toc1 R2_C1 R2_C2 R2_C3 R2_C4 R2_C5 RAAE_C1 RAAE_C2 RAAE_C3 RAAE_C4 RAAE_C5]

figure
subplot(5,1,1)
hold on
plot(C_1_out(1:80),'o');
plot(PC_p1(1:80),'r.','markersize',12);
plot(C_1_out(1:80));
hold off
legend('Calculation','Kriging','Interpreter','latex','location','NorthEastOutside','fontsize',11);
xlabel('Data points','Interpreter','latex');
ylabel('Coefficient value','interpreter','latex');
title('Fitting Accuracy of $C_1$','interpreter','latex','fontsize',14);

subplot(5,1,2)
hold on
plot(C_2_out(1:80),'o');
plot(PC_p2(1:80),'r.','markersize',12);
plot(C_2_out(1:80));
hold off
legend('Calculation','Kriging','Interpreter','latex','location','NorthEastOutside','fontsize',11);
xlabel('Data points','Interpreter','latex');
ylabel('Coefficient value','interpreter','latex');
title('Fitting Accuracy of $C_2$','interpreter','latex','fontsize',14);

subplot(5,1,3)
hold on
plot(C_3_out(1:80),'o');
plot(PC_p3(1:80),'r.','markersize',12);
plot(C_3_out(1:80));
hold off
legend('Calculation','Kriging','Interpreter','latex','location','NorthEastOutside','fontsize',11);
xlabel('Data points','Interpreter','latex');
ylabel('Coefficient value','interpreter','latex');
title('Fitting Accuracy of $C_3$','interpreter','latex','fontsize',14);

subplot(5,1,4)
hold on
plot(C_4_out(1:80),'o');
plot(PC_p4(1:80),'r.','markersize',12);
plot(C_4_out(1:80));
hold off
legend('Calculation','Kriging','Interpreter','latex','location','NorthEastOutside','fontsize',11);
xlabel('Data points','Interpreter','latex');
ylabel('Coefficient value','interpreter','latex');
title('Fitting Accuracy of $C_4$','interpreter','latex','fontsize',14);

subplot(5,1,5)
hold on
plot(C_5_out(1:80),'o');
plot(PC_p5(1:80),'r.','markersize',12);
plot(C_5_out(1:80));
hold off
legend('Calculation','Kriging','Interpreter','latex','location','NorthEastOutside','fontsize',11);
xlabel('Data points','Interpreter','latex');
ylabel('Coefficient value','interpreter','latex');
title('Fitting Accuracy of $C_5$','interpreter','latex','fontsize',14);



