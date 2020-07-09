function [C1,C2,C3,C4,C5,lambda,y_e,trace] = g(a_e,b_e,c_e,e_e,x1,x2,x3)

% Time spand
t = linspace(0,15,1500);
% Design parameters


%% Math model
% Experiment model
y_e = a_e.*x1.*cos((b_e+x2)./a_e.*t) + c_e.*log(x2).*(sin(t./x1)+sqrt(x2*t).*b_e*x3)-3+e_e;

%% Principle Component Analysis (PCA)
Yc = y_e-mean(y_e);
T = cov(Yc); % Covariance Matrix
[L eig_val] = eig(T);
eig_val = eig(T);
[lambda index_eig_val] = sort(eig_val,'descend');% lambda is the eigen value sort in order

% decreasing dimension
theta = 0;
i = 1;
while theta<0.999
    trace = sum(lambda);
    num = sum(lambda(1:i));
    theta = num/trace;
    i = i+1;
end
d = 5;

for i = 1:d
    C(:,i) = 1/sqrt(lambda(i))*Yc*L(:,index_eig_val(i));
end

lambda = lambda(1:d);

C1 = C(:,1);
C2 = C(:,2);
C3 = C(:,3);
C4 = C(:,4);
C5 = C(:,5);