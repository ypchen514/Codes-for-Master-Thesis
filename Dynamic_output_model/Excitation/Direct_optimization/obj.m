%% Objective function
function cost = obj(X)

[S_overall,var] = F_sensitivity(X);

S_1 = S_overall(1,1);
St_1 = S_overall(1,2);
S_2 = S_overall(2,1);
St_2 = S_overall(2,2);
S_3 = S_overall(3,1);
St_3 = S_overall(3,2);

% cost = ((St_1-S_1)/S_1 + (S_2+S_3)/S_1)/var;
% cost = ((St_2-S_2)/S_2 + (S_1+S_3)/S_2)/var;
 cost = ((St_3-S_3)/S_3 + (S_2+S_1)/S_3)/var;


