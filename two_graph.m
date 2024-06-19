function [GS]=two_graph(DLNS,dLNS,interaction,y_train,alpha,beta)
% DLNS ---------ҩ��������
% dLNS---- -----����������
% interaction---ҩ��-������ʼ��������A��interaction=A��ֻ�ڸ�˹���������õ��ˣ�����ͼ��ɢ��û���õ�
% y_train-------ҩ��-�������������Ѿ���WKNKN�㷨
% alpha---------ƽ�����
% beta----------ƽ�����

% ���
% GS-----���յ�Ԥ�����־���

[m,n]=size(interaction);
%% ��˹���໥����
[DS_G,dS_G] = gaussiansimilarity(interaction,m,n);
% DS_G---ҩ���ĸ�˹���ƺ�
% dS_G---������ĸ�˹���ƺ�
for i=1:m
    for j=1:m
        if DLNS(i,j)==0
            DLNS(i,j)=DS_G(i,j);
        else
            DLNS(i,j)=(DLNS(i,j)+DS_G(i,j))/2;
        end     
    end
end
for i=1:n
    for j=1:n
        if dLNS(i,j)==0
            dLNS(i,j)=dS_G(i,j);
        else
             dLNS(i,j)=(dLNS(i,j)+dS_G(i,j))/2;
        end     
    end
end
dS= dLNS;    %������ģ����Ϻ�ļ���-�����ۺ������Ծ���
DS=DLNS;     %ҩ���ģ����Ϻ��ҩ��-ҩ���ۺ������Ծ���
%%����ͼ��ɢ
C=y_train';  
A_d=dS*C;    %�����ĳ�ʼ������ɢ���ʾ���
A_D=C*DS;    %ҩ��ĳ�ʼ������ɢ���ʾ���
S_d_1=zeros(n,n);
S_d_2=zeros(n,m);
S_d1=zeros(n,n);
S_d2=zeros(n,n);
S_d3=zeros(n,m);
S_d4=zeros(n,m);
nd=n;   %����---֢����
nm=m;   %����---��ҩ��
%%��һ��
for i=1:nd         % ֢����
    for j=1:nm     % ��ҩ��
   %     S_d1(:,i)=S_d1(:,i)+(A_D(i,j)*C(:,j))/sum(A_D(:,j));
   %  �Ľ��˷�ĸ,�ô��ǿ��Ա����ĸΪ0������Դ�����޹�����ϵ�����ֿ������������Ʒ����ҩ�����Ƽ���
        S_d1(:,i)=S_d1(:,i)+(A_D(i,j)*C(:,j))*2/(sum(A_D(:,j))+sum(A_D(i,:)));

%           S_d1(:,i)=S_d1(:,i)+(A_m(i,j)*C(:,j))/(((sum(A_m(:,j)))^b)*(sum(A_m(i,:)))^(1-b));
    end
    for k=1:nm
    %    S_d2(:,i)=S_d2(:,i)+(A_d(i,k)*C(:,k))/sum(A_d(:,k));
    %  �Ľ��˷�ĸ
        S_d2(:,i)=S_d2(:,i)+(A_d(i,k)*C(:,k))*2/(sum(A_d(:,k))+sum(A_d(i,:)));
%           S_d2(:,i)=S_d2(:,i)+(A_d(i,k)*C(:,k))/(((sum(A_d(:,k)))^b)*(sum(A_d(i,:)))^(1-b));
    end
    S_d_1(:,i)=alpha*S_d1(:,i)+(1-alpha)*S_d2(:,i);
end
%%�ڶ���
for i=1:nm
    for j=1:nd
     %   S_d3(:,i)=S_d3(:,i)+(A_D(j,i)*S_d_1(:,j))/sum(A_D(j,:));
      %  �Ľ��˷�ĸ
        S_d3(:,i)=S_d3(:,i)+(A_D(j,i)*S_d_1(:,j))*2/(sum(A_D(j,:))+sum(A_D(:,i)));

%         S_d3(:,i)=S_d3(:,i)+(A_m(j,i)*S_d_1(:,j))/(((sum(A_m(j,:)))^b)*(sum(A_m(:,i)))^(1-b));
    end
    for k=1:nd
       % S_d4(:,i)=S_d4(:,i)+(A_d(k,i)*S_d_1(:,k))/sum(A_d(k,:));
        S_d4(:,i)=S_d4(:,i)+(A_d(k,i)*S_d_1(:,k))*2/(sum(A_d(k,:))+sum(A_d(:,i)));

%         S_d4(:,i)=S_d4(:,i)+(A_d(k,i)*S_d_1(:,k))/(((sum(A_d(k,:)))^b)*(sum(A_d(:,i)))^(1-b));
    end
    S_d_2(:,i)=beta*S_d3(:,i)+(1-beta)*S_d4(:,i);
end
GS=S_d_2';






