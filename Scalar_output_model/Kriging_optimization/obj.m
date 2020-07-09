%% Objective Function
function s_cost = obj(x)
global S1_krig S2_krig S3_krig S4_krig St1_krig St2_krig St3_krig St4_krig
S1 = f_predictkrige(x,S1_krig);
S2 = f_predictkrige(x,S2_krig);
S3 = f_predictkrige(x,S3_krig);
S4 = f_predictkrige(x,S4_krig);

St_1 = f_predictkrige(x,St1_krig);
St_2 = f_predictkrige(x,St2_krig);
St_3 = f_predictkrige(x,St3_krig);
St_4 = f_predictkrige(x,St4_krig);

s_cost = (St_4-S4)/S4 + (S1+S3+S2)/S4;