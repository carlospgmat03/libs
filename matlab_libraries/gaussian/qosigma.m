function[out]=qosigma(N)
% gives the sigma matrix for N modes
out=sparse(2*N,2*N);
for k=1:N,
    out(2*k-1,2*k)=1;
    out(2*k,2*k-1)=-1;
end

