%% Global Sensitivity Analysis on PC score
% The main purpose is to find the PC score of each basis function that
% compose a sample function
% Under this method, the output of a simulation will be the PC scores
% Thus, sensitivity analysis can be done by analyze the change in
% a fusion of Pc scores by a weighted sum due to different operation
% parameters
clc
clear
close all

tic
N_sample_e = 20000;
%% Model parameters
% Operation parameters
x1 = 10;
x2 = 0.5;
x3 = 10;
% Random Parameters by Low-discrepancy Sampling
Q = sobolset(8,'skip',200000);% create a different sobol sequence
Q_ss = net(Q,N_sample_e);
P_a = (Q_ss(:,1)-0.5)*0.1+1; % uncertainty of a
P_b = (Q_ss(:,2)-0.5)*0.1+1; % uncertainty of b
P_c = (Q_ss(:,3)-0.5)*0.1+1; % uncertainty of c
P_n = norminv(Q_ss(:,4),0,0.01); % sensing error e_e


Q_a = (Q_ss(:,4)-0.5)*0.1+1; % uncertainty
Q_b = (Q_ss(:,5)-0.5)*0.1+1; % uncertainty
Q_c = (Q_ss(:,6)-0.5)*0.1+1; % uncertainty
Q_n = norminv(Q_ss(:,8),0,0.01); % sensing error e_e

P = [P_a P_b P_c P_n];
Q = [Q_a Q_b Q_c Q_n];
% g calculation
[PC_p1,PC_p2,PC_p3,PC_p4,PC_p5,lambda_g,y_e] = g(P(:,1),P(:,2),P(:,3),P(:,4),x1,x2,x3);% execute p
[PC_q1,PC_q2,PC_q3,PC_q4,PC_q5,lambda] = g(Q(:,1),Q(:,2),Q(:,3),P(:,4),x1,x2,x3);% execute q
var_g1 = 1/N_sample_e*(sum(PC_p1.^2))-(1/N_sample_e*sum(PC_p1))^2;
var_g2 = 1/N_sample_e*(sum(PC_p2.^2))-(1/N_sample_e*sum(PC_p2))^2;
var_g3 = 1/N_sample_e*(sum(PC_p3.^2))-(1/N_sample_e*sum(PC_p3))^2;
var_g4 = 1/N_sample_e*(sum(PC_p4.^2))-(1/N_sample_e*sum(PC_p4))^2;
var_g5 = 1/N_sample_e*(sum(PC_p5.^2))-(1/N_sample_e*sum(PC_p5))^2;

for i = 1:4
    PP = P;
    QQ = Q;
    PP(:,i) = QQ(:,i);
    R(i) = struct('number',i,'input',PP);
    [g1_R(:,i) g2_R(:,i) g3_R(:,i) g4_R(:,i) g5_R(:,i) lambda] =g(R(i).input(:,1),R(i).input(:,2),R(i).input(:,3),R(i).input(:,4),x1,x2,x3);
    % Sensitivity of C1
    num_C1(i) = struct('no',i,'value',1/N_sample_e*(sum(PC_q1.*(g1_R(:,i)-PC_p1))));
    num_D1(i) = struct('no',i,'value',1/(2*N_sample_e)*sum((PC_p1-g1_R(:,i)).^2));
    S_C1(i) = struct('design_variable',i,'main_effect',abs(num_C1(i).value/var_g1),'total_effect',abs(num_D1(i).value/var_g1));
    % Sensitivity of C2
    num_C2(i) = struct('no',i,'value',1/N_sample_e*(sum(PC_q2.*(g2_R(:,i)-PC_p2))));
    num_D2(i) = struct('no',i,'value',1/(2*N_sample_e)*sum((PC_p2-g2_R(:,i)).^2));
    S_C2(i) = struct('design_variable',i,'main_effect',abs(num_C2(i).value/var_g2),'total_effect',abs(num_D2(i).value/var_g2));
    % Sensitivity of C3
    num_C3(i) = struct('no',i,'value',1/N_sample_e*(sum(PC_q3.*(g3_R(:,i)-PC_p3))));
    num_D3(i) = struct('no',i,'value',1/(2*N_sample_e)*sum((PC_p3-g3_R(:,i)).^2));
    S_C3(i) = struct('design_variable',i,'main_effect',abs(num_C3(i).value/var_g3),'total_effect',abs(num_D3(i).value/var_g3));
    % Sensitivity of C4
    num_C4(i) = struct('no',i,'value',1/N_sample_e*(sum(PC_q4.*(g4_R(:,i)-PC_p4))));
    num_D4(i) = struct('no',i,'value',1/(2*N_sample_e)*sum((PC_p4-g4_R(:,i)).^2));
    S_C4(i) = struct('design_variable',i,'main_effect',abs(num_C4(i).value/var_g4),'total_effect',abs(num_D4(i).value/var_g4));
    % Sensitivity of C5
    num_C5(i) = struct('no',i,'value',1/N_sample_e*(sum(PC_q5.*(g5_R(:,i)-PC_p5))));
    num_D5(i) = struct('no',i,'value',1/(2*N_sample_e)*sum((PC_p5-g5_R(:,i)).^2));
    S_C5(i) = struct('design_variable',i,'main_effect',abs(num_C5(i).value/var_g5),'total_effect',abs(num_D5(i).value/var_g5));
end
total_var = sum(lambda_g);
S_1 = (S_C1(1).main_effect*lambda_g(1)+S_C2(1).main_effect*lambda_g(2)+S_C3(1).main_effect*lambda_g(3)+S_C4(1).main_effect*lambda_g(4)+S_C5(1).main_effect*lambda_g(5))/total_var;
S_2 = (S_C1(2).main_effect*lambda_g(1)+S_C2(2).main_effect*lambda_g(2)+S_C3(2).main_effect*lambda_g(3)+S_C4(2).main_effect*lambda_g(4)+S_C5(2).main_effect*lambda_g(5))/total_var;
S_3 = (S_C1(3).main_effect*lambda_g(1)+S_C2(3).main_effect*lambda_g(2)+S_C3(3).main_effect*lambda_g(3)+S_C4(3).main_effect*lambda_g(4)+S_C5(3).main_effect*lambda_g(5))/total_var;
S_n = (S_C1(4).main_effect*lambda_g(1)+S_C2(4).main_effect*lambda_g(2)+S_C3(4).main_effect*lambda_g(3)+S_C4(4).main_effect*lambda_g(4)+S_C5(4).main_effect*lambda_g(5))/total_var;
St_1 = (S_C1(1).total_effect*lambda_g(1)+S_C2(1).total_effect*lambda_g(2)+S_C3(1).total_effect*lambda_g(3)+S_C4(1).total_effect*lambda_g(4)+S_C5(1).total_effect*lambda_g(5))/total_var;
St_2 = (S_C1(2).total_effect*lambda_g(1)+S_C2(2).total_effect*lambda_g(2)+S_C3(2).total_effect*lambda_g(3)+S_C4(2).total_effect*lambda_g(4)+S_C5(2).total_effect*lambda_g(5))/total_var;
St_3 = (S_C1(3).total_effect*lambda_g(1)+S_C2(3).total_effect*lambda_g(2)+S_C3(3).total_effect*lambda_g(3)+S_C4(3).total_effect*lambda_g(4)+S_C5(3).total_effect*lambda_g(5))/total_var;
St_n = (S_C1(4).total_effect*lambda_g(1)+S_C2(4).total_effect*lambda_g(2)+S_C3(4).total_effect*lambda_g(3)+S_C4(4).total_effect*lambda_g(4)+S_C5(4).total_effect*lambda_g(5))/total_var;


S_overall = [S_1 St_1;S_2 St_2;S_3 St_3;S_n St_n];

toc







