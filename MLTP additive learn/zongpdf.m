clear
tic
clc
%% 设置参数
a1 = -2.5;    b1 = 2.5;
x0 = -2.45;   x1 = 2.45;
h = 1/50;       
h1 = (b1-a1)/2*h;
T0 = 0;       T2 = 10;
dt = 0.5*h1^2*10;
x = x0:h1:x1;  M = length(x);
N=round((T2-T0)/dt); t=T0:dt:T2;
uu2=zeros(M,1);
x_star=zeros(N+1,1);
x_star(1)=-2;   x_star(end)=2;
d=1.0094;        % d is squre Brownian paremeter
eps=0;                  % Levy paremeter
%% 循环求解
for n=1:N-1
    u1=Finalpdf(eps,d,x_star(1),T0,t(n+1));
%     m=length(u1);
%     h1=2*(x1-x0)/m;
    uu1=u1(1:1:1+M-1);      % x=a1+N*h1
    for j=1:M
        u2=Finalpdf(eps,d,x(j),t(n+1),T2);
        uu2(j)=u2(89);      % x_star(end)=x0+h1*N
    end
    U=uu1.*uu2;
%     XU(:,n)=U;
    [ind3]=find(U==(max(U)));
    x_star(n+1)=x(min(ind3));
end
toc

save doublea25x1sta1eps01alph175T1  t x_star u1 u2 uu1 uu2 x
plot(t,x_star,'r*-')
xlabel('t');
ylabel('X_m(t)');
title 'The Maximum likelihood transition path from -2 to 2（learned） '

