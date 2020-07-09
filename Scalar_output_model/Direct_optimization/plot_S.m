%% Math Model
function [S1,S2,S3,S4,St_1,St_2,St_3,St_4] = plot_S(x)
    %% Parameters
    % fixed system parameters
    alpha = 1;
    % tunable system parameters
    p1 = x(1);
    p2 = x(2);
    p3 = x(3);
    p4 = x(4);
    
    % define design variables
    x1_lb = 0;
    x1_ub = 1;
    x1_delta = x1_ub-x1_lb;
    x2_lb = 0;
    x2_ub = 1;
    x2_delta = x2_ub-x2_lb;
    x3_lb = 0;
    x3_ub = 1;
    x3_delta = x3_ub-x3_lb;
    x4_lb = 0;
    x4_ub = 1;
    x4_delta = x4_ub-x4_lb;
    % uncertainty
    N = 1e7; % take N times sampling
    Q = sobolset(8,'skip',200000);% create a different sobol sequence
    Q_ss = net(Q,N);
    %
    %% sampling for kriging
    % Q_ss_1 = lhsdesign(N,4,'smooth','off');
    % x_1 = Q_ss_1(:,1)
    %
    % Q_ss = lhsdesign(N,8,'smooth','off');
    
    P1 = x1_lb+(Q_ss(:,1)*x1_delta);
    P2 = x2_lb+(Q_ss(:,2)*x2_delta);
    P3 = x3_lb+(Q_ss(:,3)*x3_delta);
    P4 = x4_lb+(Q_ss(:,4)*x4_delta);
    Q1 = x1_lb+(Q_ss(:,5)*x1_delta);
    Q2 = x2_lb+(Q_ss(:,6)*x2_delta);
    Q3 = x3_lb+(Q_ss(:,7)*x3_delta);
    Q4 = x4_lb+(Q_ss(:,8)*x4_delta);
    
    
    
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
%     toc
    % Call indice
    St_1 = Sj(1).total_effect;
    St_2 = Sj(2).total_effect;
    St_3 = Sj(3).total_effect;
    St_4 = Sj(4).total_effect;
    S1 = Sj(1).main_effect;
    S2 = Sj(2).main_effect;
    S3 = Sj(3).main_effect;
    S4 = Sj(4).main_effect;
    % Cost function
    s_cost = (St_4-S4)/S4 + (S1+S2+S3)/S4;
    
    end