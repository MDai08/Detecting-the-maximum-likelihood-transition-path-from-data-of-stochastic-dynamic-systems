%% �������ݼ�����Ԥ����
clear;
clc;
V = load('data.mat');
X = V.P';
Y = V.Q';
L = length(find(isnan(X)));     % ����Ԥ����
X(find(isnan(X))) = [];
Y(find(isnan(Y))) = [];
K = 1000;                         % bining���ܸ���
A = ones(K-L,1);
% H = [A X X.^2 X.^3 X.^4 X.^5 X.^6 X.^7 X.^8 X.^9 X.^10 sin(X) cos(X) sin(6*X) cos(6*X) sin(11*X) cos(11*X) tanh(10*X) -10*tanh(10*X).^2+10*exp(-50*X.^2)];
% H = [A X X.^2 X.^3];
H = [A X X.^2 X.^3 X.^4 X.^5];
% H = [A X X.^2 X.^3 X.^4 X.^5 X.^6 X.^7 X.^8];
C = pinv((H'*H))*H'*Y             % ��׼��С���˷�

%% CV+SSR
[M,J] = size(C);              % M���ܵĲ�������
N = 10;                        % k�۽�����
Delta_Vec = zeros(M,1);
E_Vec = zeros(M,N);
indices = crossvalind('Kfold', H(1:K-L,M), N);       % k�۽���
fp = fopen('par_data.txt', 'wt+');                    % ��¼�м�����
for q = 1:M                                    %q�ǲ���Ϊ0�ĸ���
    G = 0;
    for k=1:N                                                %������֤k=10��10����������Ϊ���Լ�
        test = (indices == k);                                %���test��Ԫ�������ݼ��ж�Ӧ�ı��
        train = ~test;                                        %train��Ԫ�صı��Ϊ��testԪ�صı��
        train_dataH = H(train,:);                             %�����ݼ��л��ֳ�train����������
        train_dataY = Y(train,:);
        test_dataH = H(test,:);                               %test������
        test_dataY = Y(test,:);
     % SSR
        F = pinv((train_dataH'*train_dataH))*train_dataH'*train_dataY;
        for l = 1:q  
            j = find(abs(F)==min(abs(F)));
            train_dataH(:,j) = [];
            test_dataH(:,j) = [];
            F = pinv((train_dataH'*train_dataH))*train_dataH'*train_dataY;
        end
            %��¼�м����
            [m,n] = size(F);
            for i = 1:m
                if i == m
                    fprintf(fp, '%f\r\n', F(i));
                else
                    fprintf(fp, '%f  ', F(i));
                end
            end
%         end       
     % ������֤�÷� 
        E = norm(test_dataY-test_dataH*F,2);
        E_Vec(q,k) = E;
        G = G + E; 
    end
    Delta = G/k;
    Delta_Vec(q) = Delta;
%     Delta_Vec(q) = Delta;
end
fclose(fp);
save D.mat Delta_Vec M
plot(1:M-1, flipud(Delta_Vec(1:M-1)), 'b-*');%��������
axis([0.5 5.5 0 15]);
xlabel('n')
ylabel('Delta(n)')
title 'The drift term (Cross validation score) '

% w = E_Vec(6,:);
% ww = w/sum(w);
% a = [2.465728, 7.562607, -1.941752
% 2.282454, 7.601588, -1.933104
% 2.592688, 7.459268, -1.933412
% 2.156663, 7.961053, -2.000311
% 1.756445, 7.994224, -1.964990
% 2.691473, 7.153874, -1.885609
% 2.300453, 7.562811, -1.926339
% 2.610015, 7.152156, -1.874613
% 3.108890, 6.883305, -1.866419
% 2.852275, 7.273342, -1.922243
% 
% ];
% % b(1) = ww*a(:,1);
% % b(2) = ww*a(:,2);
% % b(3) = ww*a(:,3);
% b(1) = mean(a(:,1));
% b(2) = mean(a(:,2));
% b(3) = mean(a(:,3));
% % b(4) = mean(a(:,4));
% % b(5) = mean(a(:,5));
% b