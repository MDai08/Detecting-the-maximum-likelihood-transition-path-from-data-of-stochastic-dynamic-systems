function [data_x] = derWENOr2_minus(data, dx)
%
% Calculates the derivative for left-biased stencil using
% 3rd-order accurate WENO scheme
% takes 1-D data
% data: input data
%   data(3)=f^+(X(2))*p(X(2)), X(2)=-1+h,
%   data(end-2)=f^+(X(2J))*p(X(2J)), X(2J)=1-h
% data_x(1) = data_x(X(2)) (deriv of 'data')
% dx: grid resolution
% Note: before entering this function, data need to be 
% extended by 2 at the beginning and end (values don't matter)
%
% Author: Xiaofan Li

% constant parameters
ep = 1.d-6;


% output variable
data_x = zeros(length(data)-4,1);

% data(3)=f(u(x(1)), x(1)=-1, data(end-2)=f(u(x(J))=f(u(1-h))
% extrapolate the beginning and end points of data
% data(2) = 2*data(3)-data(4);
% data(1) = 2*data(2)-data(3);
% data(end-1) = 2*data(end-2)-data(end-3);
% data(end) = 2*data(end-1)-data(end-2);

% absorbing BC
data(2)=0; data(1)=0;
data(end-1)=0; data(end)=0;

% data(2)=3*data(3)-3*data(4)+data(5); 
% data(1)=3*data(2)-3*data(3)+data(4);
% data(end-1)=3*data(end-2)-3*data(end-3)+data(end-4);
% data(end)=3*data(end-1)-3*data(end-2)+data(end-3);

%Generate the divided difference tables, D1(j)=Delta+f[j],
%            D2(j) = Delta-(Delta+ f[j]))
D1 = zeros(size(data)); D2=zeros(size(data));
D1(1:end-1) = (data(2:end)-data(1:end-1));
D2(2:end-1) = data(3:end)-2*data(2:end-1)+data(1:end-2);

for j=1:(length(data)-4)

    i = j+2;
    rm = (ep + D2(i-1)^2)/(ep + D2(i)^2);
    wm = 1/(1+2*rm^2);
    data_x(j) = ((D1(i-1)+D1(i)) - wm*(D1(i-2)-2*D1(i-1)+D1(i)))/2;

end

