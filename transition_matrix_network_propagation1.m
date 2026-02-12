%calculate the transition matrix of symmetric matrix
function [M]=transition_matrix_network_propagation1(g)
t = cputime;
l=length(g);
s=sum(g,2);
G=triu(g);
M=zeros(l,l);
%M=sparse(l,l);
for i=1:l 
    z=find(G(i,:)~=0);
    for j=1:length(z)
        M(i,z(j))=g(i,z(j))/sqrt(s(i)*s(z(j)));
        %M(z(j),i)=g(i,z(j))/sqrt(s(i)*s(z(j)));
        M(z(j),i)=M(i,z(j));
    end
end
%for i=1:length(M)
%    z=find(M(i,:)~=0);
%    M(i,z)=M(i,z)./sum(M(i,z));
%end