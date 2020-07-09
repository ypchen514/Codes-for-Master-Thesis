%% Objective Function
function cost = obj(x)

global S1_krig S2_krig S3_krig S4_krig S5_krig S6_krig St1_krig St2_krig St3_krig St4_krig St5_krig St6_krig Var_krig count
global S1_est S2_est S3_est S4_est S5_est S6_est St1_est St2_est St3_est St4_est St5_est St6_est var_est

count = count + 1;

S1_est(count) = abs(f_predictkrige(x,S1_krig));
S2_est(count) = abs(f_predictkrige(x,S2_krig));
S3_est(count) = abs(f_predictkrige(x,S3_krig));
S4_est(count) = abs(f_predictkrige(x,S4_krig));
S5_est(count) = abs(f_predictkrige(x,S5_krig));
S6_est(count) = abs(f_predictkrige(x,S6_krig));
St1_est(count) = abs(f_predictkrige(x,St1_krig));
St2_est(count) = abs(f_predictkrige(x,St2_krig));
St3_est(count) = abs(f_predictkrige(x,St3_krig));
St4_est(count) = abs(f_predictkrige(x,St4_krig));
St5_est(count) = abs(f_predictkrige(x,St5_krig));
St6_est(count) = abs(f_predictkrige(x,St6_krig));
var_est(count) = abs(f_predictkrige(x,Var_krig));

%% Cost function
% theta_1
% cost = 1/(var_est(count)^(1/3))*(abs(St1_est(count)-S1_est(count))/(S1_est(count))+(S2_est(count)+S3_est(count)+S4_est(count)+S5_est(count)+S6_est(count))/(S1_est(count)));
% theta_2
% cost = 1/(var_est(count)^(1/3))*(abs(St2_est(count)-S2_est(count))/(S2_est(count))+(S1_est(count)+S3_est(count)+S4_est(count)+S5_est(count)+S6_est(count))/(S2_est(count)));
% theta_3
% cost = 1/(var_est(count)^(1/3))*(abs(St3_est(count)-S3_est(count))/(S3_est(count))+(S2_est(count)+S1_est(count)+S4_est(count)+S5_est(count)+S6_est(count))/(S3_est(count)));
% theta_4
% cost = 1/(var_est(count)^(1/3))*(abs(St4_est(count)-S4_est(count))/(S4_est(count))+(S2_est(count)+S3_est(count)+S1_est(count)+S5_est(count)+S6_est(count))/(S4_est(count)));
% theta_5
% cost = 1/(var_est(count)^(1/3))*(abs(St5_est(count)-S5_est(count))/(S5_est(count))+(S2_est(count)+S3_est(count)+S4_est(count)+S1_est(count)+S6_est(count))/(S5_est(count)));
% theta_6
cost = 1/(var_est(count)^(1/3))*(abs(St6_est(count)-S6_est(count))/(S6_est(count))+(S2_est(count)+S3_est(count)+S4_est(count)+S5_est(count)+S1_est(count))/(S6_est(count)));
