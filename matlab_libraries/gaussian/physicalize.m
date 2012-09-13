function[out] = physicalize(gamma, mode);

if exist('RandStream')
	previous_random_stream = RandStream.getDefaultStream();
end
n=size(gamma,1)/2;
[J] = symplecticForm(n);

if (mode == '2norm')
	GammaTarget=sdpvar(2*n,2*n,'symmetric','real');
	F = set([]);
	F=F + set(GammaTarget + i*J > 0);
	sol = solvesdp(F, trace((gamma-GammaTarget)*(gamma-GammaTarget)), sdpsettings( 'verbose',0));
	out=double(GammaTarget);
elseif (mode == '1norm')
	GammaTarget=sdpvar(2*n,2*n,'symmetric','real');
	A=sdpvar(2*n,2*n,'full','real');
	B=sdpvar(2*n,2*n,'full','real');
	sdpvar t;
	%t=sdpvar(1,1,'real');
	F = set([]);
	F=F + set([A (gamma-GammaTarget); (gamma-GammaTarget) B] > 0);
	F=F + set(trace(A)+trace(B) <2*t);
	F=F + set(GammaTarget + i*J > 0);
	sol = solvesdp(F, t, sdpsettings( 'verbose',0));
	% sol = solvesdp(F, t, sdpsettings( 'verbose',0));
	out=double(GammaTarget);
else 
	fprintf('not implemented yet');
	out=NaN;
end

if exist('RandStream')
	RandStream.setDefaultStream(previous_random_stream);
end
end




