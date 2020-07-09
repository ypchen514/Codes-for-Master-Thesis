function s_cost = obj_cost(x)
global S1_cost_kriging
s_cost = f_predictkrige(x,S1_cost_kriging);