
function [kd,km] = gaussiansimilarity(interaction,nd,nm)
% kd---ҩ���ĸ�˹���ƺ�
% km---������ĸ�˹���ƺ�
% interaction---��ʼ�� ҩ��-���� ��������
%A: Binary relations between disease and miRNA, 1st column:miRNA, 2nd column:disease

%calculate gamad for Gaussian kernel calculation
% norm:����interaction��f����    f����Ϊ����A��Frobenius��������Ϊ����A����Ԫ�صľ���ֵƽ�����ܺͣ��ٿ�����
 %gamad = nd/(norm(interaction,'fro')^2);%��fro����interaction��interactionת�õĻ��ĶԽ��ߺ͵�ƽ��������sqrt(sum(diag(interaction'*interaction)))
gamad = (norm(interaction,'fro')^2)/nd;%���Դ���
%calculate Gaussian kernel for the similarity between disease: kd
C=interaction;
kd=zeros(nd,nd);
%kd=gpuArray(gd); %gpuת��
%D=C*C';     %���;�����ˣ����ױ���,����ϡ������ʾ
C1 = sparse(C);
C2 = sparse(C');
D1=C1*C2;
D=full(D1);
for i=1:nd
    for j=i:nd
        kd(i,j)=exp(-gamad*(D(i,i)+D(j,j)-2*D(i,j)));%kd����һ���������ξ���
    end
end
kd=kd+kd'-diag(diag(kd));%diag(A)���ؾ���A�����Խ����ϵ���ֵ,diag(diag(A))����ֻ���о���A�����Խ����ϵ���ֵ�ĶԽǾ�������λ��Ϊ0����ʱkd��Ϊ�������ξ���

%calculate gamam for Gaussian kernel calculation

%gamam = nm/(norm(interaction,'fro')^2);
gamam = (norm(interaction,'fro')^2)/nm;%���Դ���
%calculate Gaussian kernel for the similarity between miRNA: km
km=zeros(nm,nm);
%km=gpuArray(gm); %gpuת��
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
