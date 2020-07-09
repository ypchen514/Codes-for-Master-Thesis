%% Polynomial Chaos Expansion - 1D
% This is a basic polynomial expansion for a random variable 
% assume math function is y = tanh(x)
clear
clc
close all
% Create Legendre Function
% Variable 1
syms x_1
phi_1 = [1;x_1;1/2*(3*x_1^2-1);1/2*(5*x_1^3-3*x_1);...% 0-3
    1/8*(35*x_1^4-30*x_1^2+3);1/8*(63*x_1^5-70*x_1^3+15*x_1);... % 4-5
    1/16*(231*x_1^6-315*x_1^4+105*x_1^2-5);1/16*(429*x_1^7-693*x_1^5+315*x_1^3-35*x_1)]; % 6-7


% syms x_2
% phi_2 = [1;x_2;1/2*(3*x_2^2-1);1/2*(5*x_2^3-3*x_2);...% 0-3
%     1/8*(35*x_2^4-30*x_2^2+3);1/8*(63*x_2^5-70*x_2^3+15*x_2);... % 4-5
%     1/16*(231*x_2^6-315*x_2^4+105*x_2^2-5);1/16*(429*x_2^7-693*x_2^5+315*x_2^3-35*x_2)]; % 6-7


Qs = sobolset(1);
Q = net(Qs,100);

phi_1 = phi_1(1:6,:);
%% Parameter
% Define x
x = (Q-0.5)*2;% -1 <= x <= 1
% sub random variables in polynomial
phi = subs(phi_1',x);
% Define polynomial coefficient of x
C_x = double(inv(phi'*phi)*phi'*x);
x_re = phi_1'*C_x;



%% Output
y = tanh(x);
C_y = double(inv(phi'*phi)*phi'*y);
y_re = phi_1'*C_y;

%% Validation
for i = 1:100
Q_rand = rand(1,1); % given random variable
x_real(i) = (Q_rand-0.5)*2;
x_pol(i) = double(subs(x_re,x_real(i)));

y_real(i) = tanh(x_real(i));
y_pol(i) = double(subs(y_re,x_real(i)));
end

subplot(3,1,1)
plot(x_real)
hold on
plot(x_pol);

subplot(3,1,2)
plot(y_real);
hold on
plot(y_pol);

subplot(3,1,3)
hist(y_real-y_pol);
