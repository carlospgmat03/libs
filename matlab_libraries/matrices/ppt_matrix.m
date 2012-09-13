function[out] = ppt_matrix(ModesAB)
N=sum(ModesAB);
out = sparse(1:2*N,1:2*N,1);
for k=2:2:2*ModesAB(1)
	out(k,k)=-1;
end
