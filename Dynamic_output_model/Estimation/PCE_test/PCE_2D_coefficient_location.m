%% Polynomial Chaos Expansion 2D

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

syms x_2
phi_2 = [1;x_2;1/2*(3*x_2^2-1);1/2*(5*x_2^3-3*x_2);...% 0-3
    1/8*(35*x_2^4-30*x_2^2+3);1/8*(63*x_2^5-70*x_2^3+15*x_2);... % 4-5
    1/16*(231*x_2^6-315*x_2^4+105*x_2^2-5);1/16*(429*x_2^7-693*x_2^5+315*x_2^3-35*x_2)]; % 6-7

%% Multivariate Polynomials
% The polynomial chaos expansion term are the same for input and system
% output. So, we only need to do this process once to generate

% Multiply Term
P = 7; % truncation, P is the highest order of the polynomial
phi_term1 = [0;1;2;3;4;5;6;7];
phi_term2 = [0;1;2;3;4;5;6;7];
PHI = [];
location = [];
Pos = 0;
for i = 1:length(phi_term1)
    for j = 1:length(phi_term2)
        if phi_term1(i)+phi_term2(j) <= P % This combinition will be collected
            Mul = [i-1,j-1];
            Pos = Pos+1;
            order = sum(Mul);
            if nnz(Mul) ==1 
                location = [location;[Pos,order]];
            end
            PHI_temp = phi_1(i)*phi_2(j);
            PHI = [PHI;PHI_temp];
        end
    end
end



%% Collocation Point
Qs = sobolset(2);
Q = net(Qs,100);


%% Parameter
% By giving random variable(0-1) with n dimension, we can use polynomial
% chaos expansion to infer what system input is.

% Define x
x = (Q-0.5)*2;% -1 <= x <= 1
x1 = 0.1*x(:,1)+0.06;
x2 = 6*x(:,2);
% sub random variables in polynomial
phi = double(subs(PHI',{x_1,x_2},{x(:,1),x(:,2)}));% create collocation point
% Define polynomial coefficient of x
C_x1 = double(inv(phi'*phi)*phi'*x1);
C_x2 = double(inv(phi'*phi)*phi'*x2);
x1_re = PHI'*C_x1;
x2_re = PHI'*C_x2;

%% Output
y = sin(x1)+cos(x2)+exp(x1.*x2+10)*0.0001 ;
C_y = double(inv(phi'*phi)*phi'*y);
y_re = PHI'*C_y;

%% Variance
for i = 1:length(location)
   Var_term(i,:) = [C_x2(location(i,1)),PHI(location(i,1)),location(i,2)]; 
end    

Var_y = var(x2);
Var_pol = double(sum(Var_term(:,1).^2*2./(2*Var_term(:,3)+1))/2);


%% Validation
Q_rand = (rand(100,2)-0.5)*2; % given random variable
x1_real = 0.1*Q_rand(:,1)+0.06;
x2_real = 6*Q_rand(:,2);
for i = 1:length(Q_rand)
    x1_pol(i) = double(subs(x1_re,[x_1,x_2],[Q_rand(i,1),Q_rand(i,2)]));
    x2_pol(i) = double(subs(x2_re,[x_1,x_2],[Q_rand(i,1),Q_rand(i,2)]));
    y_pol(i) = double(subs(y_re,[x_1,x_2],[Q_rand(i,1),Q_rand(i,2)]));
end
y_real = sin(x1_real)+cos(x2_real)+exp(x1_real.*x2_real+10)*0.0001;

subplot(5,1,1)
plot(x1_real)
hold on
plot(x1_pol);
xlabel('Random Variables','Interpreter','latex');
ylabel('Value','Interpreter','latex');
title('$x_1$','Interpreter','latex','Fontsize',18);
legend('Real Value','Polynomial Chaos Expansion','Interpreter','latex');

subplot(5,1,2)
plot(x2_real)
hold on
plot(x2_pol);
xlabel('Random Variables','Interpreter','latex');
ylabel('Value','Interpreter','latex');
title('$x_2$','Interpreter','latex','Fontsize',18);
legend('Real Value','Polynomial Chaos Expansion','Interpreter','latex');

subplot(5,1,3)
plot(y_real);
hold on
plot(y_pol);
xlabel('Random Variables','Interpreter','latex');
ylabel('Value','Interpreter','latex');
title('$y$','Interpreter','latex','Fontsize',18);
legend('Real Value','Polynomial Chaos Expansion','Interpreter','latex');

subplot(5,1,4)
plot(((y_real-y_pol')./y_real)*100);
xlabel('Random Variables','Interpreter','latex');
ylabel('$100\%$','Interpreter','latex');
title('Error Percentage of $y$','Interpreter','latex','Fontsize',18);


subplot(5,1,5)
hist((abs(y_real-y_pol')./y_real*100),20);
xlabel('Error Percentage $\%$','Interpreter','latex');
ylabel('Cumulatied Number','Interpreter','latex');
title('Error Percentage Distribution','Interpreter','latex','Fontsize',18);
