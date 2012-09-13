function[out] = AttenuationChannel(state,theta,modes)

out=state;
N=size(state,1)/2;

for j1=1:size(modes,2)
	mode = modes(j1);
	Gate = eye(2*N);
	Gate(2*mode,2*mode)=cos(theta);
	Gate(2*mode -1 , 2*mode -1 )=cos(theta);
	out = Gate*out*Gate';

	Gate = zeros(2*N);
	Gate(2*mode,2*mode)=sin(theta)*sin(theta);
	Gate(2*mode -1 , 2*mode -1 )=sin(theta)*sin(theta);
	out = out + Gate;
end

return
end


