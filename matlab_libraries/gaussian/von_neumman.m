function[out] = von_neumman(Gamma);


A=i*qosigma(8)*Gamma;
d=sort(real(eig(A*A)));
d=d(1:2:end);
out=0;
for j1=1:size(d,1)
	% [j1 d(j1) h(d(j1),'bosons')]
	out = out + h(d(j1),'bosons');
end



