function[out] = how_entangled(gamma, modesAB);
% Here, 
% xe=1 is a separable state. if we lower the restriction, it is more entangled. 
% is we find a very big xe implies that the state is very very special, 
% as this poses 

% >> fprintf('Of a thermal state %f.6\nAnd of a PPT jens state %f3.6\n',how_entangled( ThermalState(0.5,8), [2 2]), how_entangled( BoundEntPPT, [2 2]) )   
% Of a thermal state 2.000000.6
% And of a PPT jens state 0.8958873.6




if exist('RandStream')
	previous_random_stream = RandStream.getDefaultStream();
end
na = modesAB(1);
nb = modesAB(2);
xe=sdpvar(1,1,'symmetric');
Gamma1=sdpvar(2*na,2*na,'symmetric');
Gamma2=sdpvar(2*nb,2*nb,'symmetric');
F = set([]);
F=F + set(gamma - [Gamma1 zeros(2*na, 2*nb); zeros(2*nb, 2*na) Gamma2]   > 0);
F=F + set(Gamma1 + i*xe*qosigma(na) > 0);
F=F + set(Gamma2 + i*xe*qosigma(nb) > 0);
sol = solvesdp(F, -xe , sdpsettings( 'verbose',0));

out = double(xe);
if exist('RandStream')
	RandStream.setDefaultStream(previous_random_stream);
end
end



