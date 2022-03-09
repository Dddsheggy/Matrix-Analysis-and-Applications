function R = traversal_CSSP(A, k)
% Input Arguments:
% A: the matrix
% k: the number of columns to choose
%
% Output Arguments:
% R: the set of column indices

[~,N]=size(A); 
col_list=nchoosek(1:N,k);
num_p=nchoosek(N,k);
er=zeros(num_p,1);
for i=1:num_p
    P=A(:,col_list(i,:));
    Ps=P*pinv(P);
    e0=A-Ps*A;
    er(i)=norm(e0,'fro');   
end
[min_er,min_id]=min(er);
R=col_list(min_id,:)';
end

