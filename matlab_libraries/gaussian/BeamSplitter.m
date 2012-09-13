function[out]= BeamSplitter(Modes, TotalNumberModes, theta,phase)

% Con theta=0 es la indentidad, que representa una refleccion de 100%
% en la forma en que lo estamos viendo 
% >> BeamSplitter([1 2], 2, 0)   
% ans =
%      1     0     0     0
%      0     1     0     0
%      0     0     1     0
%      0     0     0     1



if (nargin == 2)
	theta=pi/4;
	phase=0;
end
if (nargin == 3)
	phase=0;
end
out = OperatorOnModes(SingleModeBeamSplitter(theta,phase), Modes, TotalNumberModes);
end


