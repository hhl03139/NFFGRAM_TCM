function nearst_neighbor_matrix=calculate_neighbors(distance_matrix,neighbor_num)%% calculate the nearest K neighbors
%���ﷵ�ص�nearst_neighbor_matrix�����ǣ�01���󣩣�����1��Ӧdistance_matrix��ǰneighbor_num���������ٵ����ݣ�0��ʾ����neighbor_num���Ƿ�Χ��
  [sv si]=sort(distance_matrix,2,'ascend');%2��ʾ����distance_matrix��ÿһ���ﶼ���������У���������sv��ʾ��si��ʾdistance_matrix��sv���������к��λ��������ascend����������
  [row_num,col_num]=size(distance_matrix);%��ʱrow_num=col_num
  nearst_neighbor_matrix=zeros(row_num,col_num);%����һ����СΪrow_num*row_num�Ŀհ׾���
  index=si(:,1:neighbor_num);%����distance_matrix��sv���������к��ǰneighbor_num��λ��������������distance_matrixÿ��ǰneighbor_num�����پ����Ԫ��
  for i=1:row_num
       nearst_neighbor_matrix(i,index(i,:))=1;%������������λ�÷��ص�nearst_neighbor_matrix�ϣ�������λ����1��û�б�����λ����0
  end
end
