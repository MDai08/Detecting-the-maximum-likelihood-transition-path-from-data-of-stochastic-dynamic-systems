%% 调用数据及数据预处理
clear;
clc;
W = load('doublea25x1sta1eps01alph175T3.mat');
V = load('odoublea25x1sta1eps01alph175T3.mat');

T = V.t';
X = V.x_star';
Y = W.x_star';

plot(T,X,'b.-','MarkerSize',8)
hold on
plot(T,Y,'ro-','MarkerSize',4)
legend('original','learned')
xlabel('t');
ylabel('X_m(t)');
title 'The Maximum likelihood transition path from -2 to 2 '
hold off