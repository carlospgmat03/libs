function[out] = trace_norm(A)
out = sum( sqrt( eig(A*transpose(conj(A)) )) ) ;
end

