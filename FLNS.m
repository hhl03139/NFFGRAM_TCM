function W=FLNS(feature_matrix,tag,neighbor_num,regulation) %% Using the method of label propagation to predict the interaction


%function W=Label_Propagation(feature_matrix,tag,neighbor_num,regulation) %% Using the method of label propagation to predict the interaction
    distance_matrix=calculate_instances_FDim(feature_matrix);   %�÷���ά��

    nearst_neighbor_matrix=calculate_neighbors(distance_matrix,neighbor_num);%���ﷵ�ص�nearst_neighbor_matrix�����ǣ�01���󣩣�����1��Ӧdistance_matrix��ǰneighbor_num���������ٵ����ݣ�0��ʾ����neighbor_num���Ƿ�Χ��
    W=optimization_similairty_matrix(feature_matrix,nearst_neighbor_matrix,tag,regulation);
end

