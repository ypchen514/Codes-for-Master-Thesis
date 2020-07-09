%% Kriging Optimization of Time-dependent output math model
% ===================================================
%% Fit Kriging of operation parameters
total_sample = 125;
X_ss = sobolset(3,'skip',40000);
x = net(X_ss,total_sample);
X1 = (x(:,1)-0.5)*10+5;
X2 = x(:,2)*0.5;
X3 = (x(:,3)-0.5)*10+5;
X = [X1 X2 X3];


tic
parfor i = 1:length(X)
    [S_overall,var] = F_sensitivity(X(i,:));
    S_1(i) = S_overall(1,1);
    St_1(i) = S_overall(1,2);
    S_2(i) = S_overall(2,1);
    St_2(i) = S_overall(2,2);
    S_3(i) = S_overall(3,1);
    St_3(i) = S_overall(3,2);
    Var(i) = var;
end

lb = [0,0,0];
ub = [1,1,1];
fittype = 1;
SCFtype = 1;
max_variance = 1;
run_num = 0;

% S1_krig = f_variogram_fit(x, S_1', lb, ub);
% S2_krig = f_variogram_fit(x, S_2', lb, ub);
% S3_krig = f_variogram_fit(x, S_3', lb, ub);
% St1_krig = f_variogram_fit(x, St_1', lb, ub);
% St2_krig = f_variogram_fit(x, St_2', lb, ub);
% St3_krig = f_variogram_fit(x, St_3', lb, ub);
% Var_krig = f_variogram_fit(x, Var', lb, ub);

S1_krig = f_variogram_fit(X, S_1', lb, ub);
S2_krig = f_variogram_fit(X, S_2', lb, ub);
S3_krig = f_variogram_fit(X, S_3', lb, ub);
St1_krig = f_variogram_fit(X, St_1', lb, ub);
St2_krig = f_variogram_fit(X, St_2', lb, ub);
St3_krig = f_variogram_fit(X, St_3', lb, ub);
Var_krig = f_variogram_fit(X, Var', lb, ub);

save('S1_krig.mat','S1_krig');
save('S2_krig.mat','S2_krig');
save('S3_krig.mat','S3_krig');
save('St1_krig.mat','St1_krig');
save('St2_krig.mat','St2_krig');
save('St3_krig.mat','St3_krig');
save('Var_krig.mat','Var_krig');
