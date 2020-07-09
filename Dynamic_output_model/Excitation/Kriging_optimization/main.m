%% Direct Optimization Algorithm - - Main
% main file
clc
clear

load('S1_krig.mat');
load('S2_krig.mat');
load('S3_krig.mat');
load('St1_krig.mat');
load('St2_krig.mat');
load('St3_krig.mat');
load('Var.mat');

count = 0;

global S1_krig S2_krig S3_krig St1_krig St2_krig St3_krig Var_krig count
global S1_est S2_est S3_est St1_est St2_est St3_est var_est

S_1 = S1_est;
S_2 = S2_est;
S_3 = S3_est;
St_1 = St1_est;
St_2 = St2_est;
St_3 = St3_est;
Var = var_est;


tic
% setting
options.conTol = 1e-6;
options.display = 0;
fileInfo.fName = 'obj';
% algorithm
result = UMDIRECT(fileInfo,[5,0,0],[10,0.5,10], options);
% result
x_best = result.xBest;
y_best = result.fBest;
toc

history = result.rect.x;

