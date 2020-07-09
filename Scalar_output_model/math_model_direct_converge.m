clear
clc
close all
tic
%% Parameters
% fixed system parameters
alpha = 1;
% tunable system parameters
p1 = -2;
p2 = -0.5;
p3 = -2;
p4 = -0.5;

% Sampling
N_samp = [1e3,5e3,1e4,5e4,1e5,5e5,1e5,5e6,1e7]; % take N times sampling

for samp_n = 1:length(N_samp)
    N = N_samp(samp_n);
    Q = sobolset(8);% create a different sobol sequence
    Q_ss = net(Q,N);%Sobol set for 8 dimension
    P1 = Q_ss(:,1);
    P2 = Q_ss(:,2);
    P3 = Q_ss(:,3);
    P4 = Q_ss(:,4);
    Q1 = Q_ss(:,5);
    Q2 = Q_ss(:,6);
    Q3 = Q_ss(:,7);
    Q4 = Q_ss(:,8);
    
    for j = 1:1
        %% Calculation
        % Q1 = rand(N,1);
        % Q2 = rand(N,1);
        % Q3 = rand(N,1);
        % Q4 = rand(N,1);
        
        
        P = [P1, P2, P3, P4];
        Q = [Q1, Q2, Q3, Q4];
        
        % g calculation
        g_P = g(alpha,p1,p2,p3,p4,P); % sub P into g to calculate g(P)
        g_Q = g(alpha,p1,p2,p3,p4,Q); % sub Q into g to calculate q(P)
        var_g = 1/N*(sum(g_P.^2))-(1/N*sum(g_P))^2; %Total variance
        
        for i = 1:4
            PP = P;
            QQ = Q;
            PP(:,i) = QQ(:,i);
            R(i) = struct('number',i,'in',PP); % will generate R matrix as sample matrix
            g_R(i) = struct('no',i,'value',g(alpha,p1,p2,p3,p4,R(i).in));% calulate each g(R^j)
            C(i) = struct('no',i,'value',1/N*(sum(g_Q.*(g_R(i).value-g_P))));
            D(i) = struct('no',i,'value',1/(2*N)*sum((g_P-g_R(i).value).^2));
            Sj(i) = struct('design_variable',i,'main_effect',C(i).value/var_g,'total_effect',D(i).value/var_g);
            
        end
        
        Sj_s(j) = struct('no',j,'S1',Sj(1).main_effect,'S2',Sj(2).main_effect,'S3',Sj(3).main_effect,'S4',Sj(4).main_effect);
    end
    toc
    Sj_all(samp_n) = struct('samples',N,'Sj',Sj);
    
end