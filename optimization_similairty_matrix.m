%%'regulation1':LN similarity, 'regulation2': RLN similarity
function W=optimization_similairty_matrix(feature_matrix,nearst_neighbor_matrix,tag,regulation) %%quadratic programming
   row_num=size(feature_matrix,1);%1���ؾ���feature_matrix������
   W=zeros(1,row_num);%WΪ1��row_num�е�0����
   if tag==1%
       row_num=1;
   end
   for i=1:row_num
       nearst_neighbors=feature_matrix(logical(nearst_neighbor_matrix(i,:)'),:);%��feature_matrix�����i��������ŷ����þ���ǰneighbor_num �����������  
       neighbors_num=size(nearst_neighbors,1);%����nearst_neighbors�������������û���⣬�������neighbor_num����֮ǰɸѡ�������������
       G1=repmat(feature_matrix(i,:),neighbors_num,1)-nearst_neighbors;%������feature_matrix(i,:)����neighbors_num*1�飬��G1��neighbors_num*1������feature_matrix(i,:)ƽ�̶���
       %������˼�ǽ�feature_matrix��i���������������Ӧnearst_neighbor�����������ɵľ���ͬά�ȵľ��󣬲���������֮��Ĳ�ࡣ
       G2=repmat(feature_matrix(i,:),neighbors_num,1)'-nearst_neighbors';%ͬ��
       %G2��ʵ����G1��ת�ã���������������
       G11=sparse(G1);
       G22=sparse(G2);
       if regulation=='regulation2'
      %   G_i=G1*G2+eye(neighbors_num);%eye(N),����N*N��С�ĵ�λ���󣻶Խ���Ϊ1������Ϊ0
         %�����ǽ������i������������������������Ĳ�������ڻ��ľ���regulation1�Ĳ������ڶԽ����Ƿ��1��G_i��neighbors_num*neighbors_num
          G_i=G11*G22+eye(neighbors_num);  %����ϡ�������ˣ���ʡ�ڴ�
       end
       if regulation=='regulation1'
       %  G_i=G1*G2;%������
         G_i=G11*G22;%������
       end
       H=2*G_i;%����2��
       f=[];
       A=[];
       if isempty(H)%���H����Ϊ�㣬�ͷ���A
           A;
       end
       
       b=[];
       Aeq=ones(neighbors_num,1)';%����һ��ȫ1����1*neighbors_num
       beq=1;
       lb=zeros(neighbors_num,1);%����һ��ȫ1����neighbors_num*1
       ub=[];
       options=optimset('Display','off');
       %����Ķ��ι滮������matlab������ģ�������
       [w,fval]= quadprog(H,f,A,b,Aeq,beq,lb,ub,[],options);%quadprog���ι滮
       w=w';
       W(i,logical(nearst_neighbor_matrix(i,:)))=w;%�����ι滮��������Ȩֵ������ڵ�λ���������и�ֵ    
   end
 %  H
end