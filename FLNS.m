function W=FLNS(feature_matrix,tag,neighbor_num,regulation) %% Using the method of label propagation to predict the interaction


%function W=Label_Propagation(feature_matrix,tag,neighbor_num,regulation) %% Using the method of label propagation to predict the interaction
    distance_matrix=calculate_instances_FDim(feature_matrix);   %用分形维数

    nearst_neighbor_matrix=calculate_neighbors(distance_matrix,neighbor_num);%这里返回的nearst_neighbor_matrix矩阵是（01矩阵），其中1对应distance_matrix中前neighbor_num个距离最少的数据，0表示不在neighbor_num考虑范围内
    W=optimization_similairty_matrix(feature_matrix,nearst_neighbor_matrix,tag,regulation);
end

