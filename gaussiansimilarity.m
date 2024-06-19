
function [kd,km] = gaussiansimilarity(interaction,nd,nm)
% kd---药物间的高斯相似核
% km---疾病间的高斯相似核
% interaction---初始的 药物-疾病 关联矩阵
%A: Binary relations between disease and miRNA, 1st column:miRNA, 2nd column:disease

%calculate gamad for Gaussian kernel calculation
% norm:返回interaction的f范数    f范数为矩阵A的Frobenius范数定义为矩阵A各项元素的绝对值平方的总和，再开根号
 %gamad = nd/(norm(interaction,'fro')^2);%‘fro’：interaction和interaction转置的积的对角线和的平方根，即sqrt(sum(diag(interaction'*interaction)))
gamad = (norm(interaction,'fro')^2)/nd;%测试代码
%calculate Gaussian kernel for the similarity between disease: kd
C=interaction;
kd=zeros(nd,nd);
%kd=gpuArray(gd); %gpu转换
%D=C*C';     %大型矩阵相乘，容易崩溃,改用稀疏矩阵表示
C1 = sparse(C);
C2 = sparse(C');
D1=C1*C2;
D=full(D1);
for i=1:nd
    for j=i:nd
        kd(i,j)=exp(-gamad*(D(i,i)+D(j,j)-2*D(i,j)));%kd最后成一个上三角形矩阵
    end
end
kd=kd+kd'-diag(diag(kd));%diag(A)返回矩阵A的正对角线上的数值,diag(diag(A))返回只含有矩阵A的正对角线上的数值的对角矩阵，其他位置为0，此时kd仍为上三角形矩阵

%calculate gamam for Gaussian kernel calculation

%gamam = nm/(norm(interaction,'fro')^2);
gamam = (norm(interaction,'fro')^2)/nm;%测试代码
%calculate Gaussian kernel for the similarity between miRNA: km
km=zeros(nm,nm);
%km=gpuArray(gm); %gpu转换
%E=C'*C;
E1=C2*C1;
E=full(E1);
for i=1:nm
    for j=i:nm
        km(i,j)=exp(-gamam*(E(i,i)+E(j,j)-2*E(i,j)));
    end
end
km=km+km'-diag(diag(km));
end
