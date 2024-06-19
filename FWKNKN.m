function [MD_mat_new] = FWKNKN( MD_mat, MM_mat, DD_mat, K, r )
%输入
   % MD_mat---药物-疾病的关联矩阵
   % MM_mat---药物-药物相似性矩阵
   % DD_mat---疾病-疾病相似性矩阵
   % K--K近邻的K，邻域的大小
   % r--衰减因子
% 输出
   % MD_mat_new---更新后的药物-疾病关联矩阵



[rows,cols]=size(MD_mat);
y_m=zeros(rows,cols);  
y_d=zeros(rows,cols);  

knn_network_m = KNN( MM_mat, K );  %for miRNA  对药物
for i = 1 : rows   %药物的个数
         w=zeros(1,K);
        [sort_m,idx_m]=sort(knn_network_m(i,:),2,'descend'); 
        sum_m=sum(sort_m(1,1:K));   
        for j = 1 : K
            w(1,j)=r^(j-1)*sort_m(1,j); 
            y_m(i,:) =  y_m(i,:)+ w(1,j)* MD_mat(idx_m(1,j),:); 
        end                      
            y_m(i,:)=y_m(i,:)/sum_m;              
end

knn_network_d = KNN( DD_mat , K );  %for disease  对疾病
for i = 1 : cols      %疾病的个数
        w=zeros(1,K);
        [sort_d,idx_d]=sort(knn_network_d(i,:),2,'descend');
        sum_d=sum(sort_d(1,1:K));
        for j = 1 : K
            w(1,j)=r^(j-1)*sort_d(1,j);
            y_d(:,i) =  y_d(:,i)+ w(1,j)* MD_mat(:,idx_d(1,j)); 
        end                      
            y_d(:,i)=y_d(:,i)/sum_d;               
end

a1=1;
a2=1;
y_md=(y_m*a1+y_d*a2)/(a1+a2);  

for i = 1 : rows
    for j = 1 : cols
        MD_mat_new(i,j)=max(MD_mat(i,j),y_md(i,j));   %更新的药物疾病关联矩阵

    end
end

end



