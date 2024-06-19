function distance_matrix=calculate_instances_FDim(feature_matrix)
%matrix2Tensor

%计算图片的分形维数
cellmax = 1024;
[row_num,col_num]=size(feature_matrix);%feature_matrix矩阵的维度row_num*col_num
distance_matrix=zeros(row_num,row_num);%构建一个大小为row_num*row_num的空白矩阵
D=zeros(row_num,1);
for i=1:row_num  %feature_matrix的行数
    D(i) = FractalDim(feature_matrix(i,:),cellmax);   %与上面的效果相同，速度快
end
for i = 1:row_num
    for j=i+1:row_num
        distance_matrix(i,j) = abs(D(i)-D(j));
    end
end
% 让距离矩阵为对称矩阵
distance_matrix = triu(distance_matrix,0) + tril(distance_matrix',-1);
%将自己置最大值，表示最近邻不考虑自身向量
distance_matrix = distance_matrix + eye(row_num,row_num)*10*max(max(distance_matrix));

%以下程序用于对每一行i的距离进行升序排序，找到每个药物（或疾病）的最相关的前K个药物
[distance_matrix_sort,I] = sort(distance_matrix,2);
end