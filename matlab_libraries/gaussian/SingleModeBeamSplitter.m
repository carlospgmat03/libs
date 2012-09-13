function[out]= SingleModeBeamSplitter(theta,phase)
% A single mode phase gate
if (nargin == 0)
	theta=pi/4;
end
% out = [1 0 1 0;0 1 0 1;-1 0 1 0; 0 -1 0 1]/sqrt(2);
out = [cos(theta) 0 -sin(theta) 0;...
	0 cos(theta) 0 -sin(theta);...
	sin(theta) 0 cos(theta) 0;...
	0 sin(theta) 0 cos(theta)];


if (nargin>1)
	out=out*WhichModePhaseGate(phase, 1, 2);
	out=out*WhichModePhaseGate(phase, 2, 2);	
end

end

