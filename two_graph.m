function [GS]=two_graph(DLNS,dLNS,interaction,y_train,alpha,beta)
% DLNS ---------药物相似性
% dLNS---- -----疾病相似性
% interaction---药物-疾病初始关联矩阵A，interaction=A，只在高斯核相似性用到了，二分图扩散中没有用到
% y_train-------药物-疾病关联矩阵，已经过WKNKN算法
% alpha---------平衡参数
% beta----------平衡参数

% 输出
% GS-----最终的预测评分矩阵

[m,n]=size(interaction);
%% 高斯核相互作用
[DS_G,dS_G] = gaussiansimilarity(interaction,m,n);
% DS_G---药物间的高斯相似核
% dS_G---疾病间的高斯相似核
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
dS= dLNS;    %疾病间的，整合后的疾病-疾病综合相似性矩阵
DS=DLNS;     %药物间的，整合后的药物-药物综合相似性矩阵
%%二部图扩散
C=y_train';  
A_d=dS*C;    %疾病的初始特征扩散概率矩阵
A_D=C*DS;    %药物的初始特征扩散概率矩阵
S_d_1=zeros(n,n);
S_d_2=zeros(n,m);
S_d1=zeros(n,n);
S_d2=zeros(n,n);
S_d3=zeros(n,m);
S_d4=zeros(n,m);
nd=n;   %列数---症候数
nm=m;   %行数---草药数
%%第一步
for i=1:nd         % 症候数
    for j=1:nm     % 草药数
   %     S_d1(:,i)=S_d1(:,i)+(A_D(i,j)*C(:,j))/sum(A_D(:,j));
   %  改进了分母,好处是可以避免分母为0（无资源，即无关联关系），又可以提高冷门物品（草药）的推荐率
        S_d1(:,i)=S_d1(:,i)+(A_D(i,j)*C(:,j))*2/(sum(A_D(:,j))+sum(A_D(i,:)));

%           S_d1(:,i)=S_d1(:,i)+(A_m(i,j)*C(:,j))/(((sum(A_m(:,j)))^b)*(sum(A_m(i,:)))^(1-b));
    end
    for k=1:nm
    %    S_d2(:,i)=S_d2(:,i)+(A_d(i,k)*C(:,k))/sum(A_d(:,k));
    %  改进了分母
        S_d2(:,i)=S_d2(:,i)+(A_d(i,k)*C(:,k))*2/(sum(A_d(:,k))+sum(A_d(i,:)));
%           S_d2(:,i)=S_d2(:,i)+(A_d(i,k)*C(:,k))/(((sum(A_d(:,k)))^b)*(sum(A_d(i,:)))^(1-b));
    end
    S_d_1(:,i)=alpha*S_d1(:,i)+(1-alpha)*S_d2(:,i);
end
%%第二步
for i=1:nm
    for j=1:nd
     %   S_d3(:,i)=S_d3(:,i)+(A_D(j,i)*S_d_1(:,j))/sum(A_D(j,:));
      %  改进了分母
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






