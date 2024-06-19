function [auc_apur_10,pre_rec,formula_pre] = HerbRecommendation_attentionLOSS_KDHR(K,et,kD,kd,alpha,beta,topK)
% [auc_apur_10,pre_rec,formula_pre] = HerbRecommendation_attentionLOSS_KDHR(9,0.4,20,50,0.1,0.6,10)
% K,et -----------FWKNKN( A, DCS, dSS, K, et);
% kD   -----------DLNS=FLNS(y_train,0,kD,'regulation2'); 药物间的线性邻域相似性
% kd--------------dLNS=FLNS(y_train',0,kd,'regulation2'); 证候间的线性邻域相似性
% alpha，beta-----two_graph(DLNS,dLNS,A,y_train,alpha,beta)

%KDHR805_390数据集
dSS = xlsread('KDHR/ss_onehot_matrix.csv');  
A = xlsread('KDHR/sh_390_onehot_matrix.csv');  
A=A';                                       
DCS  = xlsread('KDHR/hh_onehot_matrix.csv');
%去掉标签行和列
dSS(1,:) = [];
dSS(:,1) = [];
DCS(1,:) = [];
DCS(:,1) = [];
A(1,:) = [];
A(:,1) = [];
A_ori=A;%805*390

y_train=FWKNKN( A, DCS, dSS, K, et); 
dLNS=FLNS(y_train',0,kd,'regulation2');    
DLNS=FLNS(y_train,0,kD,'regulation2');    

%计算稀疏率
% 初始药物-疾病关联矩阵的稀疏率
RateSpar_A = length(find(A==0))/(length(A(:,1))*length(A(1,:))); %矩阵中0的个数除以元素总数
%经FWKNKN运算后的稀疏率
RateSpar_Awk =length(find(y_train==0))/(length(y_train(:,1))*length(y_train(1,:)));
F_1_ori=two_graph(DLNS,dLNS,A,y_train,alpha,beta);  %F_1_ori 评分矩阵
%经FWKNKN运算,再经two_graph后的稀疏率
RateSpar_Awk_two =length(find(F_1_ori<=0.01))/(length(F_1_ori(:,1))*length(F_1_ori(1,:)));
pre_label_score = F_1_ori(:);
label_y = A_ori(:);
i=1;
auc = roc_1(pre_label_score,label_y,'red');
aupr =pr_cure(pre_label_score,label_y,'blue');

F_1_ori_ori=F_1_ori;
index=find(A_ori==1);
auc = zeros(1,10);
aupr=zeros(1,10);
%十折交叉验证
for i = 1:10
    
    indices = crossvalind('Kfold', length(index), 10);
    A = A_ori;
    F_1_ori=F_1_ori_ori;
    for cv = 1:10
        cv;
        index_2 = find(cv == indices);
        A(index(index_2)) = 0;          
        y_train_1=FWKNKN( A, DCS, dSS, K, et );
        dLNS1=FLNS(y_train_1',0,kd,'regulation2');
        DLNS1=FLNS(y_train_1,0,kD,'regulation2');
        F_1=two_graph( DLNS1,dLNS1,A,y_train_1,alpha,beta);
        F_1_ori(index(index_2)) = F_1(index(index_2));
        A = A_ori;
    end
    pre_label_score = F_1_ori(:);
    label_y = A_ori(:);
    auc(i) = roc_1(pre_label_score,label_y,'red');
    aupr(i) =pr_cure(pre_label_score,label_y,'blue');
    auc(i)
end
auc_ave = mean(auc);
auc_std = std(auc);
aupr_ave = mean(aupr);
aupr_std = std(aupr);

auc_apur_10 = [auc_ave,auc;aupr_ave,aupr];
auc_apur_10 = auc_apur_10'

RateSpar_Awk_two_tenfold =length(find(F_1_ori<=0.01))/(length(F_1_ori(:,1))*length(F_1_ori(1,:)));


test_users = xlsread('KDHR\pres_data\pres_data/pS_onehot_test.csv');  % 药方-证候0-1矩阵
test_users(1,:)=[];
test_users(:,1)=[];
test_herb=xlsread('KDHR\PH_PS\PH_PS/PH_testset.csv');    %原药方中的草药配伍,来自 test.txt
test_herb(1,:)=[];
test_herb(:,1)=[];

Num_herb = sum(~isnan(test_herb'));  %每个药方含有的草药数


%用矩阵的稀疏率
RateSpar_old = RateSpar_Awk_two_tenfold;
RateSpar_new = RateSpar_old+0.1;
%用 （原始得分矩阵-新得分矩阵）的范数作为损失函数
Norm_old = norm(A_ori-F_1_ori,'fro');
Norm_new = norm(A_ori-F_1_ori,'fro')+10;


LOSS = abs(RateSpar_new-RateSpar_old)*abs(Norm_new-Norm_old); %稀疏率差值*范数差值，作为损失函数

while (LOSS > 0.001)
    Q=F_1_ori;
    K=Q;
    V=Q;
    score=Q*K'/sqrt(length(Q(1,:)));  %缩放点积模型
    score_exp = exp(score);             %矩阵中每个元素求e的次方
    score_exp_sum = sum(score_exp,2);   %每行求和
    SoftMax=score_exp./score_exp_sum;
    Score_pred = SoftMax*V;
    F_1_ori = Score_pred;

    RateSpar_old = RateSpar_new;
    RateSpar_new =length(find(F_1_ori<=0.01))/(length(F_1_ori(:,1))*length(F_1_ori(1,:)));

    Norm_old = Norm_new;
    Norm_new = norm(A_ori-F_1_ori,'fro');
    LOSS = abs(RateSpar_new-RateSpar_old)*abs(Norm_new-Norm_old);
end


as_sh = xlsread('KDHR\KDHR_data\KDHR_data/as_sh.csv');
as_sh(1,:) = [];
as_sh(:,1) = [];

%草药出现频数
idx_fre_herb = xlsread('KDHR\KDHR_data\KDHR_data/idx_fre_herb.csv');
fre_as_sh = xlsread('KDHR\KDHR_data\KDHR_data/fre_as_sh.csv');
fre_as_sh(1,:) = [];
fre_as_sh(:,1) = [];

as_sh1 = fre_as_sh./(repmat(idx_fre_herb(:,3)',length(as_sh(:,1)),1));

%对药物重要性指标as_sh1进行softmax
score_exp = exp(as_sh1);             %矩阵中每个元素求e的次方
score_exp_sum = sum(score_exp,2);   %每行求和
as_sh1_SoftMax=score_exp./score_exp_sum;

S_H_new = F_1_ori'.*as_sh.*as_sh1_SoftMax;  % 点乘 比 加 效果好
formula_pre_score=test_users*S_H_new;
[score_sort,formula_pre] = sort(formula_pre_score,2,"descend");

%matlab得到的草药编号-1=原test中草药编号
formula_pre=formula_pre-ones(size(formula_pre));
%再求查准率和查全率
%topK=5;
for i=1:length(formula_pre)
    formula_K{i} = intersect(formula_pre(i,1:topK),test_herb(i,:));
    precision_K(i) = length(formula_K{i})/topK;
    recall_K(i)= length(formula_K{i})/Num_herb(i);
end
precision_K_ave = mean(precision_K)
recall_K_ave = mean(recall_K)
F1_score_K= 2*precision_K_ave* recall_K_ave/(precision_K_ave+ recall_K_ave)


pre_rec = [precision_K_ave,precision_K;recall_K_ave,recall_K];  %查准率和查全率，第一行为平均值，从第二行开始为每个药方的查准率和查全率
pre_rec = pre_rec';

end
