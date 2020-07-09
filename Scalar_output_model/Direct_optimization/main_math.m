%% main
clc
clear


tic
% setting
options.conTol = 1e-3;
options.display = 0;
fileInfo.fName = 'obj';
% algorithm
result = UMDIRECT(fileInfo,[-5,-5,-5,-5],[0,0,0,0], options);
% result
x_best = result.xBest;
y_best = result.fBest;
toc