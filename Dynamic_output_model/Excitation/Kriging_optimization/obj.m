%% Objective function
function cost = obj(x)

global S1_krig S2_krig S3_krig St1_krig St2_krig St3_krig Var_krig count
global S1_est S2_est S3_est St1_est St2_est St3_est var_est

count = count + 1;

S1_est(count) = f_predictkrige(x,S1_krig);
S2_est(count) = f_predictkrige(x,S2_krig);
S3_est(count) = f_predictkrige(x,S3_krig);
St1_est(count) = f_predictkrige(x,St1_krig);
St2_est(count) = f_predictkrige(x,St2_krig);
St3_est(count) = f_predictkrige(x,St3_krig);
var_est(count) = f_predictkrige(x,Var_krig);

cost = ((St1_est(count)-S1_est(count))/S1_est(count) + (S2_est(count)+S3_est(count))/S1_est(count))/var_est(count);
% cost = ((St2_est(count)-S2_est(count))/S2_est(count) + (S1_est(count)+S3_est(count))/S2_est(count))/var_est(count);
%  cost = ((St3_est(count)-S3_est(count))/S3_est(count) + (S2_est(count)+S1_est(count))/S3_est(count))/var_est(count);
 
 
