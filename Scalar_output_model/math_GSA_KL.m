%% Global Sensitivity Analysis on PC score
% The main purpose is to find the PC score of each basis function that
% compose a sample function
% Under this method, the output of a simulation will be the PC scores
% Thus, sensitivity analysis can be done by analyze the change in
% a fusion of Pc scores by a weighted sum due to different operation
% parameters

clear
close all


%% Observe Target
Ms = 1;
p = 0;
mu = 0;

N_sample_e = 4000;

%% Model parameters
% Random Parameters by Low-discrepancy Sampling
Q = sobolset(3,'skip',2000000);% create a different sobol sequence
Q_ss = net(Q,N_sample_e);
a_e = (Q_ss(:,1)-0.5)*0.1+1; % uncertainty
b_e = (Q_ss(:,2)-0.5)*0.1+1; % uncertainty
c_e = (Q_ss(:,3)-0.5)*0.1+1;
e_e = normrnd(0,0.01,N_sample_e,1); % experimental error








