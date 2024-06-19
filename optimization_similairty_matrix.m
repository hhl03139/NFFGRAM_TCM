%%'regulation1':LN similarity, 'regulation2': RLN similarity
function W=optimization_similairty_matrix(feature_matrix,nearst_neighbor_matrix,tag,regulation) %%quadratic programming
   row_num=size(feature_matrix,1);%1返回矩阵feature_matrix的行数
   W=zeros(1,row_num);%W为1行row_num列的0矩阵
   if tag==1%
       row_num=1;
   end
   for i=1:row_num
       nearst_neighbors=feature_matrix(logical(nearst_neighbor_matrix(i,:)'),:);%将feature_matrix中与第i行向量的欧几里得距离前neighbor_num 个行向量抽出  
       neighbors_num=size(nearst_neighbors,1);%返回nearst_neighbors矩阵的行数，如没意外，这里等于neighbor_num，即之前筛选的最近向量个数
       G1=repmat(feature_matrix(i,:),neighbors_num,1)-nearst_neighbors;%将向量feature_matrix(i,:)复制neighbors_num*1块，即G1由neighbors_num*1块向量feature_matrix(i,:)平铺而成
       %这里意思是将feature_matrix第i行向量构建与其对应nearst_neighbor个最近向量组成的矩阵同维度的矩阵，并计算他们之间的差距。
       G2=repmat(feature_matrix(i,:),neighbors_num,1)'-nearst_neighbors';%同理
       %G2其实就是G1的转置，方便计算下面操作
       G11=sparse(G1);
       G22=sparse(G2);
       if regulation=='regulation2'
      %   G_i=G1*G2+eye(neighbors_num);%eye(N),返回N*N大小的单位矩阵；对角线为1，其他为0
         %这里是将上面的i行向量矩阵与其最近邻向量的差异矩阵内积的矩阵，regulation1的差异在于对角线是否加1，G_i是neighbors_num*neighbors_num
          G_i=G11*G22+eye(neighbors_num);  %利用稀疏矩阵相乘，节省内存
       end
       if regulation=='regulation1'
       %  G_i=G1*G2;%看上面
         G_i=G11*G22;%看上面
       end
       H=2*G_i;%乘以2，
       f=[];
       A=[];
       if isempty(H)%如果H矩阵为零，就返回A
           A;
       end
       
       b=[];
       Aeq=ones(neighbors_num,1)';%构建一个全1矩阵，1*neighbors_num
       beq=1;
       lb=zeros(neighbors_num,1);%构建一个全1矩阵，neighbors_num*1
       ub=[];
       options=optimset('Display','off');
       %下面的二次规划函数，matlab会崩溃的，？？？
       [w,fval]= quadprog(H,f,A,b,Aeq,beq,lb,ub,[],options);%quadprog二次规划
       w=w';
       W(i,logical(nearst_neighbor_matrix(i,:)))=w;%将二次规划构建出的权值按最近邻的位置索引进行赋值    
   end
 %  H
end