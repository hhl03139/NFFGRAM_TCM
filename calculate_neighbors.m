function nearst_neighbor_matrix=calculate_neighbors(distance_matrix,neighbor_num)%% calculate the nearest K neighbors
%这里返回的nearst_neighbor_matrix矩阵是（01矩阵），其中1对应distance_matrix中前neighbor_num个距离最少的数据，0表示不在neighbor_num考虑范围内
  [sv si]=sort(distance_matrix,2,'ascend');%2表示矩阵distance_matrix的每一行里都按升序排列，排序结果用sv表示，si表示distance_matrix按sv排序结果排列后的位置索引，ascend是升序排序
  [row_num,col_num]=size(distance_matrix);%此时row_num=col_num
  nearst_neighbor_matrix=zeros(row_num,col_num);%构建一个大小为row_num*row_num的空白矩阵
  index=si(:,1:neighbor_num);%保留distance_matrix按sv排序结果排列后的前neighbor_num个位置索引，即保留distance_matrix每行前neighbor_num个最少距离的元素
  for i=1:row_num
       nearst_neighbor_matrix(i,index(i,:))=1;%将保留索引的位置返回到nearst_neighbor_matrix上，保留的位置是1，没有保留的位置是0
  end
end
