%药物-药物间的分形关联性，
feature_matrix = y_train;    % 草药-证候的初始关联矩阵
distance_matrix=calculate_instances_FDim(feature_matrix);
%中，对distance_matrix的每一行从高到低排序，得到相应编号，就是与药物i关联最近的前k个药物，得到药物-药物关联关系
%类似的，可得证候-证候关联关系
%与利用公式从原始药方中挖掘得到的关联关系可以做一对比


%以下程序用于对每一行i的距离进行升序排序，找到每个药物（或疾病）的最相关的前K个药物
[relationshipScore_fractal_sort,relationship_fractal] = sort(distance_matrix,2);

%欧式距离
distance_matrix=calculate_instances(feature_matrix);
[relationshipScore_Euclid_sort,relationship_Euclid] = sort(distance_matrix,2);

%相关系数
distance_matrix=calculate_instances_corr(feature_matrix);
[relationshipScore_corr_sort,relationship_corr] = sort(distance_matrix,2);


