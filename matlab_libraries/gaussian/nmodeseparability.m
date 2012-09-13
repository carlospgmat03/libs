function[out] = nmodeseparability(Gamma);
% We check if the Gamma in is in all its modes separable
if exist('RandStream')
	previous_random_stream = RandStream.getDefaultStream();
end
ntotal=length(Gamma)/2;
F = set([]);
Gsep=cell(ntotal,1);
for j=1:ntotal
  	Gsep(j)={sdpvar(2,2,'symmetric')};
end
GammaSeparable =  DirectSum(Gsep);
F=F + set(Gamma - GammaSeparable   > 0);
xe=sdpvar(1,1,'symmetric');
F=F + set(GammaSeparable + i*xe*qosigma(ntotal) > 0);
sol = solvesdp(F, -xe , sdpsettings( 'verbose',0));
out = double(xe);
if exist('RandStream')
	RandStream.setDefaultStream(previous_random_stream);
end
