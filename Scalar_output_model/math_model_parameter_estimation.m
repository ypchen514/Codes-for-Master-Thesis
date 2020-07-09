%% Parameter Estimation of math model (Scalar output)
clc
clear
% Simulate the real output value under optimal operation parameter
% to excite p_1, p_2, use [-1.1366,-0.1857,-4.7374,-4.0052]
% to excite p_3, p_4, use [-0.0944,-0.6667,-0.0833,-0.6833]

%% uncertain parameter
% uncertainty
N = 1e7; % take N times sampling
Q = sobolset(4,'skip',200000);% create a different sobol sequence
Q_ss = net(Q,N);

p1 = Q_ss(:,1);
p2 = Q_ss(:,2);
p3 = Q_ss(:,3);
p4 = Q_ss(:,4);


global x_12 x_34 y_real_12 y_real_34 x_1234 y_real_1234 x_all y_real_all
% p = [0.5234,0.1331,0.1621,0.054];
p = rand(4,1);
x_validate = [-1.5,-0.8,-0.23,-2];
x_12 = [-1.3272,-0.0309,-4.9074,-4.9691];
x_34 = [-0.0926,-0.6481,-0.0926,-0.6481];
x_1234 = [-2,-0.5,-2,-0.5];
x_all = [-5,-3,-2,-1];

y_real_12 = g(1,p(1),p(2),p(3),p(4),x_12);
y_real_34 = g(1,p(1),p(2),p(3),p(4),x_34);
y_real_1234 = g(1,p(1),p(2),p(3),p(4),x_1234);
y_real_all = g(1,p(1),p(2),p(3),p(4),x_all);


%% round 1
% Only adjust p3,p4
tic
% setting
options.conTol = 1e-10;
options.display = 0;
options.maxIter = 10000;
% options.maxfCount= 500;
fileInfo.fName = 'obj';
% algorithm
result = UMDIRECT(fileInfo,[0,0,0,0],[1,1,1,1], options);
% result
p_best1 = result.xBest;
y_best1 = result.fBest;
toc

% %% round 2
% % let p3,p4 moves only a little, focus on adjusting p1,p2
% p_ub = [min(1,p_best1(1)+0.5),min(1,p_best1(2)+0.5),min(1,p_best1(3)+0.1),min(1,p_best1(4)+0.1)];
% p_lb = [max(0,p_best1(1)-0.5),max(0,p_best1(2)-0.5),max(0,p_best1(3)-0.1),max(0,p_best1(4)-0.1)];
% 
% options.conTol = 1e-7;
% options.display = 0;
% options.maxIter = 10000;
% % options.maxfCount=500;
% fileInfo.fName = 'obj2';
% % algorithm
% result = UMDIRECT(fileInfo,[0,0,p_lb(3:4)],[1,1,p_ub(3:4)], options);
% % result
% p_best2 = result.xBest;
% y_best2 = result.fBest;
% toc
% 
% 
% %% round 3
% p_ub = min(p_best2+0.2,1);
% p_lb = max(p_best2-0.2,0);
% 
% options.conTol = 1e-7;
% options.display = 0;
% options.maxIter = 10000;
% % options.maxfCount=500;
% fileInfo.fName = 'obj3';
% % algorithm
% result = UMDIRECT(fileInfo,p_lb,p_ub, options);
% % result
% p_best3 = result.xBest;
% y_best3 = result.fBest;
% toc

% est_g = g(1,p_best3(1),p_best3(2),p_best3(3),p_best3(4),x_validate);
% real_g = g(1,p(1),p(2),p(3),p(4),x_validate);

est_g = g(1,p_best1(1),p_best1(2),p_best1(3),p_best1(4),x_validate);
real_g = g(1,p(1),p(2),p(3),p(4),x_validate);

