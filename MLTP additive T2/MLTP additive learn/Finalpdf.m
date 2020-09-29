%%%ï¿½ï¿½T0ï¿½ï¿½ï¿½ãµ½T1
function U=Finalpdf(eps,d,x0,T1,T2)
 alpha=1.75;
 a1=-2.5;b1=2.5;
 eps =( eps^alpha)*(2/(b1-a1))^alpha;
 d = d^2*4/(b1-a1)^2;
%  f = @(x,t)       (k*(0.5*s+0.2*s*tanh((x-265)/10)-eplon*sig*x.^4))*2/(b1-a1);
 f= @(x,t)        (-1.0225*x.^3+4.0362*x)*2/(b1-a1);
 rr=1;             %right
 rl=1;             %left
 
 h=1/50;
 J=rr/h;             %right
 L=rl/h;             %left
 dt = 0.5*h^2*10;
 dtdx = dt/h;
%  X=-rl:h:rr;
 XX = a1:(b1-a1)/2*h:b1;
 nft=round((T2-T1)/dt);
 t=T1:dt:T1+nft*dt;
 Jt=L+J-1;      %Jt=total number of unknowns  
 ff=zeros(nft+1,1);
 for j=1:nft+1
     ff(j)=max(abs(f(XX,t(j))));
 end
 fpmax=max(ff);
 x=-(rl+rr):h:(rl+rr);
 %ÉèÖÃ²ÎÊý
 C=-zeta(alpha-1)*h^(2-alpha);          % correction term u''(x)
 cons=alpha*gamma((1+alpha)/2)/(2^(1-alpha)*sqrt(pi)*gamma(1-alpha/2));  %the constant of alpha-stable symmetric jump measure
%  coeff for the diffusion and the correction term
 Chh=( d/2 + cons*C*eps )/h^2;
 c1=eps*cons/alpha;
 c2=eps*cons*h;
 
 
 
 b=zeros(Jt,1);
 a=zeros(Jt,1);
 c=zeros(Jt,1); 
 
 
 %nonintegral part nature
 % coefficient of U_j
 b(2:Jt-1) =( -2*Chh- c1*(1./(x(J+3:2*J+L-1)+rl).^alpha+1./(rr-x(J+3:2*J+L-1)).^alpha))';
 % one-sided diff near boundaries
 b(1)  = 2*Chh- c1*(1/(rl+x(J+2))^alpha+1/(rr-x(J+2))^alpha); 
 b(Jt) = 2*Chh- c1*(1/(rl+x(2*J+L))^alpha+1/(rr-x(2*J+L))^alpha); % one-sided diff
 % %nonintegral part/absorbing
 % % coefficient of U_j
 % b(2:Jt-1) = (-2*Chh- c1*(1./(x(J+3:2*J+L-1)+rl).^alpha+1./(rr-x(J+3:2*J+L-1)).^alpha))';
 % % one-sided diff near boundaries
 % b(1)  =2*Chh- c1*(1/(rl+x(J+2))^alpha+1/(rr-x(J+2))^alpha); 
 % b(Jt) = 2*Chh- c1*(1/(rl+x(2*J+L))^alpha+1/(rr-x(2*J+L))^alpha); % one-sided diff
 
 
 a= Chh*ones(Jt,1);  % coefficient of U_(j-1)
 c= Chh*ones(Jt,1);  % coefficient of U_(j+1) 
 c(1)  = -5*Chh; % one-sided diff
 a(Jt) = -5*Chh; % one-sided diff
 vp2 = zeros(Jt,1); vp2(3) = 4*Chh;  % one-sided diff
 vp3 = zeros(Jt,1); vp3(4) =  -Chh; % one-sided diff 
 vm2 = zeros(Jt,1); vm2(Jt-2) = 4*Chh;  % one-sided diff
 vm3 = zeros(Jt,1); vm3(Jt-3) =  -Chh; % one-sided diff 
 % integral part
 for j=-L+1:J-1
     b(j+L)= b(j+L) - c2*( sum(1./abs(x(J+2-j:L+J)).^(1+alpha)) ...
                       + sum(1./abs(x(L+J+2:2*J+L-j)).^(1+alpha)) ...
         + .5/abs(x(J+1-j))^(1+alpha) + .5/abs(x(2*J+L+1-j))^(1+alpha) );  
 end
 A=spdiags([vm3 vm2 [a(2:end); 0] b ...
           [0; c(1:end-1)] vp2 vp3],-3:3,Jt,Jt);
 % coefficient of u_(j+k) 
 j=-L+1;
 br1=c2*[1./abs(x(J+2-j:L+J)).^(1+alpha)  0  1./abs(x(L+J+2:2*J+L-j)).^(1+alpha)];
 B=toeplitzmultaux(br1',br1);
%B=zeros(size(A));
%for j=-L+1:J-1
%  B(L+j,:)=[1./abs(x(J+2-j:L+J)).^(1+alpha)  0  1./abs(x(L+J+2:2*J+L-j)).^(1+alpha)];
%end
%   XX=-rl:h:rr;
XX = a1 : (b1-a1)/2*h :b1;
UU=sqrt(40/pi)*exp(-40*(XX-x0).^2); %gaussian

% UU= T0./(pi.*(X.^2+T0^2));     %cauchy
%UU=ones(size(X))./(2*r);
U=UU(2:end-1)';
Un=U; U1=U; U2=U;
%U=UU(2:end-1);
nu = length(U);
data = zeros(nu+4,1);
x_star=zeros(nft+1,1);
x_star(1)=x0;
for nt=1:nft
    %U1 = U + dt*(A+c2*B)*U;
    U1 = U + dt*A*U + dt*toeplitzmult2(B,U);
    % global Lax-Friedrichs(LF) flux splitting
    data(3:nu+2) = (f(XX(2:end-1),t(nt)) + fpmax)'.*U/2;
    fx1 = derWENOr2_minus(data,h);
    data(3:nu+2) = (f(XX(2:end-1),t(nt)) - fpmax)'.*U/2;
    fx2 = derWENOr2_plus(data,h);
    U1 = U1 - dtdx*(fx1+fx2);
     
    %U2 = 0.75*U + U1/4 + (dt/4)*(A+c2*B)*U1;
    U2 = 0.75*U + U1/4 + (dt/4)*A*U1 + (dt/4)*toeplitzmult2(B,U1);
    data(3:nu+2) = (f(XX(2:end-1),t(nt)) + fpmax)'.*U1/2;
    fx1 = derWENOr2_minus(data,h);
    data(3:nu+2) = (f(XX(2:end-1),t(nt)) - fpmax)'.*U1/2;
    fx2 = derWENOr2_plus(data,h);
    U2 = U2 - (dtdx/4)*(fx1+fx2);
    
    %Un = U/3 + 2*U2/3 + (2*dt/3)*(A+c2*B)*U2;
    Un = U/3 + 2*U2/3 + (2*dt/3)*A*U2 + (2*dt/3)*toeplitzmult2(B,U2);
    data(3:nu+2) = (f(XX(2:end-1),t(nt)) + fpmax)'.*U2/2;
    fx1 = derWENOr2_minus(data,h);
    data(3:nu+2) = (f(XX(2:end-1),t(nt)) - fpmax)'.*U2/2;
    fx2 = derWENOr2_plus(data,h);
    Un = Un - (2*dtdx/3)*(fx1+fx2);     
    U=Un;
%         [ind3]=find(U==(max(U)));
%     x_star(nt+1)=XX(min(ind3)+1);
end
%  plot(XX(2:end-1)',U,'--')
  

 
