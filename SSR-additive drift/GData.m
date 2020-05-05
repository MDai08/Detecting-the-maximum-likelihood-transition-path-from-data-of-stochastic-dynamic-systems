clc;clear;
tic
randn('state',100);
rand('state',100);
%% Parameters
h = 0.0001;  %迭代步长
n1 = 101;  %迭代次数
n2 = 1000000; %初值个数 1e6 sample for additive noise, 3e6 sample for multiplicative noise
% sigma = ones(n1-1,n2);  %%variation

%% Input
% X = 4*rand(1,n2) - 2;  %Uniform distribution

% % adapt grid
% X = zeros(1,n2);
% X(1:2*n2/5) = -1.5:(1/(2*n2/5-1)):-0.5;
% X(2*n2/5+1:3*n2/5) = -0.5:(1/(n2/5-1)):0.5;
% X(3*n2/5+1:end) = 0.5:(1/(2*n2/5-1)):1.5;

% % adapt grid
% X = zeros(1,n2);
% X(1:2*n2/5) = -0.75:(0.5/(2*n2/5-1)):-0.25;
% X(2*n2/5+1:3*n2/5) = -0.25:(0.5/(n2/5-1)):0.25;
% X(3*n2/5+1:end) = 0.25:(0.5/(2*n2/5-1)):0.75;

X = -2:(4/(n2-1)):2;  %grid point

Z = zeros(n1,n2);
Z(1,:) = X;  %%initial
W = sqrt(h)*randn(n1-1,n2);
for i = 1:n2
    for j = 1:n1-1
        Z(j+1,i) = Z(j,i) + (4*Z(j,i) - Z(j,i)^3)*h + 1*W(j,i);  %% Euler method for additive noise
%         Z(j+1,i) = Z(j,i) + (4*Z(j,i) - Z(j,i)^3)*h + (Z(j,i)+1)*W(j,i);  %% Euler method for multiplicative noise
    end
end
Y = Z(end,:);
X1 = (Y-X)./0.01;
% X1 = (Y-X).^2./0.01;
%% binning
K = 1000;            % bining的总个数
[Z1,E] = discretize(X,K);
for k = 1:K
%     W(k) = length(find(Z==k))/length(X);
    P(k) = sum(X(find(Z1==k)))/length(X(find(Z1==k)));
    Q(k) = sum(X1(find(Z1==k)))/length(X1(find(Z1==k)));
end
% WW = diag(W);
save data.mat X Y Z E P Q %WW
