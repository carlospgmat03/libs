function[out] = ApplyBeamSplitter(Gamma, modes, AttenuationTheta, theta, phase);

if (nargin <4 )
	theta=pi/4;
end
if (nargin <5 )
	phase=0;
end

NumberOfModes = size(Gamma,2)/2;
Gate = BeamSplitter(modes, NumberOfModes,theta,phase);
out= Gate*Gamma*Gate';
if (nargin > 2)
	out = AttenuationChannel(out, AttenuationTheta, modes);
end

end
