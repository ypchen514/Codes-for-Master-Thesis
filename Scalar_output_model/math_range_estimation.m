%% Possible Range Pridiction
clear
clc
tic
% Real value of p
p_real = rand(1,4);
%% Original range of unknown parameters
N = 4000; % take N times sampling
Q = sobolset(4,'skip',200000);% create a different sobol sequence
Q_ss = net(Q,N);
p1 = Q_ss(:,1);
p2 = Q_ss(:,2);
p3 = Q_ss(:,3);
p4 = Q_ss(:,4);

%% Situation 0
x_1234 = [-2,-0.5,-2,-0.5];
% Real value
y_real_0 = g(1,p_real(1),p_real(2),p_real(3),p_real(4),x_1234);

% Fit Kriging
for i = 1:N
    y_x0(i) = g(1,p1(i),p2(i),p3(i),p4(i),x_1234);
end

lb = [0,0,0,0];
ub = [1,1,1,1];
fittype = 1;
SCFtype = 1;

input = Q_ss;
output = y_x0';
kparam_0 = f_variogram_fit(input, output, lb, ub);
toc_fit = toc

% Scan Kriging
N_scan = 50000; % take N times sampling
Q = sobolset(4,'skip',200000);% create a different sobol sequence
Scan_0 = net(Q,N_scan);
y_krig_0 = f_predictkrige(Scan_0, kparam_0);
data_krig_0 = [Scan_0 abs(y_krig_0-y_real_0)];
data_sort_krig_0 = sortrows(data_krig_0,5);
data_sort_krig_0 = data_sort_krig_0(1:1000,:);
toc_predict = toc

%% Situation 1 -- excite p1,p2
x_12 = [-1.3272,-0.0309,-4.9074,-4.9691];

% Real value
y_real_1 = g(1,p_real(1),p_real(2),p_real(3),p_real(4),x_12);

% Fit Kriging
for i = 1:N
    y_x1(i) = g(1,p1(i),p2(i),p3(i),p4(i),x_12);
end

lb = [0,0,0,0];
ub = [1,1,1,1];
fittype = 1;
SCFtype = 1;

input = Q_ss;
output = y_x1';
kparam_1 = f_variogram_fit(input, output, lb, ub);
toc_fit = toc

% Scan Kriging
N_scan = 50000; % take N times sampling
Q = sobolset(4,'skip',20000000);% create a different sobol sequence
Scan_1 = net(Q,N_scan);
y_krig_1 = f_predictkrige(Scan_1, kparam_1);
data_krig_1 = [Scan_1 abs(y_krig_1-y_real_1)];
data_sort_krig_1 = sortrows(data_krig_1,5);
data_sort_krig_1 = data_sort_krig_1(1:1000,:);
toc_predict = toc

%% Situation 2 -- excite p3,p4
x_34 = [-0.0926,-0.6481,-0.0926,-0.6481];

% Real value
y_real_2 = g(1,p_real(1),p_real(2),p_real(3),p_real(4),x_34);

% Fit Kriging
for i = 1:N
    y_x2(i) = g(1,p1(i),p2(i),p3(i),p4(i),x_34);
end

lb = [0,0,0,0];
ub = [1,1,1,1];
fittype = 1;
SCFtype = 1;

input = Q_ss;
output = y_x2';
kparam_2 = f_variogram_fit(input, output, lb, ub);
toc_fit = toc

% Scan Kriging
N_scan = 50000; % take N times sampling
Q = sobolset(4,'skip',900000);% create a different sobol sequence
Scan_2 = net(Q,N_scan);
y_krig_2 = f_predictkrige(Scan_2, kparam_2);
data_krig_2 = [Scan_1 abs(y_krig_2-y_real_2)];
data_sort_krig_2 = sortrows(data_krig_2,5);
data_sort_krig_2 = data_sort_krig_2(1:1000,:);
toc_predict = toc

%% Overall
data_fusion = [Scan_0 abs((y_krig_0-y_real_0)/y_real_0)+abs((y_krig_1-y_real_1)/y_real_1)+abs((y_krig_2-y_real_2)/y_real_2)];
data_sort_fusion = sortrows(data_fusion,5);
data_sort_fusion = data_sort_fusion(1:100,:);

% subplot(2,2,1)
% hist(data_sort_fusion(:,1),20);
% subplot(2,2,2)
% hist(data_sort_fusion(:,2),20);
% subplot(2,2,3)
% hist(data_sort_fusion(:,3),20);
% subplot(2,2,4)
% hist(data_sort_fusion(:,4),20);

subplot(2,2,1)
hist(data_sort_krig_2(:,1),20);
subplot(2,2,2)
hist(data_sort_krig_2(:,2),20);
subplot(2,2,3)
hist(data_sort_krig_2(:,3),20);
subplot(2,2,4)
hist(data_sort_krig_2(:,4),20);