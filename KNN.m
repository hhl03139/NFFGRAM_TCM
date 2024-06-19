function [ knn_network ] = KNN( network , k )
    [rows, cols] = size( network );
    network= network-diag(diag(network));               %将主对角线元素减掉
    knn_network = zeros(rows, cols);
    [sort_network,idx]=sort(network,2,'descend');       %按行降序排列
    for i = 1 : rows
        knn_network(i,idx(i,1:k))=sort_network(i,1:k);  
    end
end