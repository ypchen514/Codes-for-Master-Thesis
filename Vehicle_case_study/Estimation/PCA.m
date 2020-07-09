function [C,trace,lambda,L,index_eig_val] = PCA(sim_data)
y_e = sim_data;
%% Principle Component Analysis (PCA)
Yc = y_e-mean(y_e);
T = cov(Yc); % Covariance Matrix
[L eig_val] = eig(T);
eig_val = eig(T);
[lambda index_eig_val] = sort(eig_val,'descend');% lambda is the eigen value sort in order

% decreasing dimension
theta = 0;
i = 1;
while theta<0.99
    trace = sum(lambda);
    num = sum(lambda(1:i));
    theta = num/trace;
    i = i+1;
end
d = i;

for i = 1:d
    C(:,i) = 1/sqrt(lambda(i))*Yc*L(:,index_eig_val(i));
end

lambda = lambda(1:d);

% C1 = C(:,1);
% C2 = C(:,2);
% C3 = C(:,3);
% C4 = C(:,4);
% C5 = C(:,5);
