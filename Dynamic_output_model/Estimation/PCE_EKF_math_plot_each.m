%% Polynomial Chaos Expansion Based Kalman -- Time Dependent Output Case
% We use a time-dependent output model as an example 
% y_e = a_e.*x1.*cos((b_e+x2)./a_e.*t) + c_e.*log(x2).*(sin(t./x1)+sqrt(x2*t).*b_e*x3)-3+e_e;
clear
clc
close all
% Create Legendre Function
% Variable 1
syms x_1
phi_1 = [1;x_1;1/2*(3*x_1^2-1);1/2*(5*x_1^3-3*x_1);...% 0-3
    1/8*(35*x_1^4-30*x_1^2+3);1/8*(63*x_1^5-70*x_1^3+15*x_1);... % 4-5
    1/16*(231*x_1^6-315*x_1^4+105*x_1^2-5);1/16*(429*x_1^7-693*x_1^5+315*x_1^3-35*x_1)]; % 6-7
subs(phi_1,0.5);

% Variable 2
syms x_2
phi_2 = [1;x_2;1/2*(3*x_2^2-1);1/2*(5*x_2^3-3*x_2);...% 0-3
    1/8*(35*x_2^4-30*x_2^2+3);1/8*(63*x_2^5-70*x_2^3+15*x_2);... % 4-5
    1/16*(231*x_2^6-315*x_2^4+105*x_2^2-5);1/16*(429*x_2^7-693*x_2^5+315*x_2^3-35*x_2)]; % 6-7

% Variable 3
syms x_3
syms x_3
phi_3 = [1;x_3;1/2*(3*x_3^2-1);1/2*(5*x_3^3-3*x_3);...% 0-3
    1/8*(35*x_3^4-30*x_3^2+3);1/8*(63*x_3^5-70*x_3^3+15*x_3);... % 4-5
    1/16*(231*x_3^6-315*x_3^4+105*x_3^2-5);1/16*(429*x_3^7-693*x_3^5+315*x_3^3-35*x_3)]; % 6-7



%% Multivariate Polynomials
% The polynomial chaos expansion term are the same for input and system
% output. So, we only need to do this process once to generate

% Multiply Term
P = 7; % truncation, P is the highest order of the polynomial
phi_term1 = [0;1;2;3;4;5;6;7];
phi_term2 = [0;1;2;3;4;5;6;7];
phi_term3 = [0;1;2;3;4;5;6;7];
PHI = [];
location = [];
Pos = 0;
for i = 1:length(phi_term1)
    for j = 1:length(phi_term2)
        for l = 1:length(phi_term3)
        if phi_term1(i)+phi_term2(j)+phi_term3(l) <= P
            Mul = [i-1,j-1,l-1];
            Pos = Pos+1;
            order = sum(Mul);
            if nnz(Mul) ==1 
                location = [location;[Pos,order]];
            end
            PHI_temp = phi_1(i)*phi_2(j)*phi_3(l);
            PHI = [PHI;PHI_temp]; % final multivariate polynomial equation 
        end
        end
    end
end

%% Collocation Point
Qs = sobolset(3);
Q = net(Qs,3*length(PHI)); % 0~1


%% Coefficient
% By giving random variable(0-1) with n dimension, we can use polynomial
% chaos expansion to infer what system input is.

% Define x
x = (Q-0.5)*2; % -1 <= x <= 1, this is the random variable for PCE estimation
x1 = ((x(:,1)/2+0.5)-0.5)*0.1+1;
x2 = ((x(:,2)/2+0.5)-0.5)*0.1+1;
x3 = ((x(:,3)/2+0.5)-0.5)*0.1+1;

% sub Collocation point in polynomial 
phi = double(subs(PHI',{x_1,x_2,x_3},{x(:,1),x(:,2),x(:,3)})); % create collocation point

% Define polynomial coefficient of x by Least Square
C_x1 = double(inv(phi'*phi)*phi'*x1);
C_x2 = double(inv(phi'*phi)*phi'*x2);
C_x3 = double(inv(phi'*phi)*phi'*x3);

% Reconstruct model input in PCE form adding calculated coefficient
x1_re = PHI'*C_x1;
x2_re = PHI'*C_x2;
x3_re = PHI'*C_x3;

%% Output
% Generate system output with uncertain model parameters
p1 = 9.8333;
p2 = 0.4883;
p3 = 0.3333;
time_space = 1500;
terminal_time = 15;
t = linspace(0,terminal_time,time_space);
y_e = x1.*p1.*cos((x2+p2)./x1.*t) + x3.*log(p2).*(sin(t./p1)+sqrt(p2*t).*x2*p3)-3;
figure
plot(t,y_e);


% Coefficient of polynomial in y_e(t_k)
% Calculate the coefficient of y_e(t_k), k = 0~k
% desire_piece = 100;
% k = terminal_time/desire_piece;
% for i = 1:(desire_piece+1)
%     t_k = k*(i-1);
%     y_e = x1.*p1.*cos((x2+p2)./x1.*t_k) + x3.*log(p2).*(sin(t_k./p1)+sqrt(p2*t_k).*x2*p3)-3+normrnd(0,0.01,length(Q),1);
%     C_y = double(inv(phi'*phi)*phi'*y_e);
%     C_y_k(i).k = i;
%     C_y_k(i).C_y = C_y;
%     y_re(i).pol = PHI'*C_y;
% end



%% Validation
% err_per = [];
%     Q_rand = rand(100,3); % random variables, 0~1
%     x = (Q_rand-0.5)*2; % -1 <= x <= 1, this is the random variable for PCE estimation
%     x1 = ((x(:,1)/2+0.5)-0.5)*0.1+1;
%     x2 = ((x(:,2)/2+0.5)-0.5)*0.1+1;
%     x3 = ((x(:,3)/2+0.5)-0.5)*0.1+1;
% for i = 1:(desire_piece+1)
%     t_k = k*(i-1);
%     y_real = x1.*p1.*cos((x2+p2)./x1.*t_k) + x3.*log(p2).*(sin(t_k./p1)+sqrt(p2*t_k).*x2*p3)-3;
%     for j = 1:length(Q_rand)
%         y_pol(j) = double(subs(y_re(i).pol,[x_1,x_2,x_3],[x(j,1),x(j,2),x(j,3)]));
%     end
%     figure(i)
%     plot(y_real);
%     hold on
%     plot(y_pol);
%     err = abs((y_real-y_pol')./y_real)*100;
%     err_per = [err_per;err];
% end
% figure
% hist(err_per,50);

%% Extended Kalman Filter
desire_piece = 20;
k = terminal_time/desire_piece;
% Real model
% assume x1_real = 0.95, x2_real = 1.01, x3_real = 1.04
x1_real = 1.015;
x2_real = 1.03;
x3_real = 0.96;

y_exp = x1_real.*p1.*cos((x2_real+p2)./x1_real.*t) + x3_real.*log(p2).*(sin(t./p1)+...
        sqrt(p2*t).*x2_real*p3)-3+normrnd(0,0.1,1,time_space);
    
% Kalman Filter
% 1. Simulation under the predict parameters
% 2. Predict covariance matrix
% 3. Estimate coefficient of each parameters

H = 1;
R = 0.01;
for j = 1:desire_piece+1
   % Simulation to desire time
   t_k = k*(j-1);
%    for n = 1:length(x)
%        x1_kf(n) = double(subs(x1_re,[x_1,x_2,x_3],[x(n,1),x(n,2),x(n,3)]));
%        x2_kf(n) = double(subs(x2_re,[x_1,x_2,x_3],[x(n,1),x(n,2),x(n,3)]));
%        x3_kf(n) = double(subs(x3_re,[x_1,x_2,x_3],[x(n,1),x(n,2),x(n,3)]));
%    end 
   x1_kf = phi*C_x1;
   x2_kf = phi*C_x2;
   x3_kf = phi*C_x3;
  
   y_sim_kf = x1_kf.*p1.*cos((x2_kf+p2)./x1_kf.*t_k) + x3_kf.*log(p2).*(sin(t_k./p1)+...
              sqrt(p2*t_k).*x2_kf*p3)-3+normrnd(0,0.01,length(Q),1);
    
   % Real model output at t = t_k
   z_k = interp1(t,y_exp,t_k);
   
   % PCE of y_sim_kf
    C_y_kf = double(inv(phi'*phi)*phi'*y_sim_kf);
   
   % Calculate covariance matrix
   P = cov([y_sim_kf,x1_kf,x2_kf,x3_kf]);
   Pyy = P(1,1);
   P_x1_y = P(1,2);
   P_x2_y = P(1,3);
   P_x3_y = P(1,4);
   
   % Update coefficient of each polynomial at a time
   for i = 1:length(PHI)
       % Kernal Delta function
       if i ==1
           delta = 1;
       else 
           delta = 0;
       end
       
       C_x1_a(i) = C_x1(i)+P_x1_y*H'*inv(R+H*Pyy*H)*(z_k*delta-H*C_y_kf(i));
       C_x2_a(i) = C_x2(i)+P_x2_y*H'*inv(R+H*Pyy*H)*(z_k*delta-H*C_y_kf(i));
       C_x3_a(i) = C_x3(i)+P_x3_y*H'*inv(R+H*Pyy*H)*(z_k*delta-H*C_y_kf(i));
   end
       C_x1 = C_x1_a';
       C_x2 = C_x2_a';
       C_x3 = C_x3_a';
       
       
       x1_final = phi*C_x1;
       x2_final = phi*C_x2;
       x3_final = phi*C_x3;
       
       x1_save(j).result = x1_final;
       x2_save(j).result = x2_final;
       x3_save(j).result = x3_final;
end    

% Final Value
x1_final = phi*C_x1;
x2_final = phi*C_x2;
x3_final = phi*C_x3;

x1_target = x1_real*ones(0.2*length(x1_final),1);
x2_target = x2_real*ones(0.2*length(x2_final),1);
x3_target = x3_real*ones(0.2*length(x3_final),1);


%% Plot result
subplot(1,3,1)
hist([x1,x1_final],100);
hold on
plot([x1_real x1_real],[0,200]);
mean(x1_final);
xlabel('$x_1$','interpreter','latex','fontsize',16)
legend('Initial Distribution','Estimate by KF','Real Value','Interpreter','latex')


subplot(1,3,2)
hist([x2,x2_final],100);
hold on
plot([x2_real x2_real],[0,200]);
mean(x2_final);
xlabel('$x_2$','interpreter','latex','fontsize',16)
legend('Initial Distribution','Estimate by KF','Real Value','Interpreter','latex')


subplot(1,3,3)
hist([x3,x3_final],100);
hold on
plot([x3_real x3_real],[0,200]);
mean(x3_final); 
xlabel('$x_3$','interpreter','latex','fontsize',16)
legend('Initial Distribution','Estimate by KF','Real Value','Interpreter','latex')

suptitle('Parameter Estimation with Extended Kalman Filter')

figure
y_val = mean(x1_final).*p1.*cos((mean(x2_final)+p2)./mean(x1_final).*t) + mean(x3_final).*log(p2).*(sin(t./p1)+...
              sqrt(p2*t).*mean(x2_final)*p3)-3;
plot(y_exp);
hold on
plot(y_val);

% figure
% for i = 1:j
%     hold on
%     histf(x3_save(i).result,20);
% end
% legend('1','2','3','4','5','6','7','8','9','10','11');
% plot([x3_real x3_real],[0,200]);
