function[out]= InitialStateProduct(Squeezing, Noise)
% These two are in general vectors, that contain the information
% 


NModes=size(Squeezing,2);
if ( NModes ~= size(Noise,2)) 
	error('initialStateproduct:ArgumentsNotEqualSize', ...
		'Arguments must have the same: %i, %i.', ...
		size(Squeezing,2), size(Noise,2))
end

out = [];
for k=1:NModes
	gamma_local = Squeezed(Squeezing(k)) + Noise(k)*eye(2);
	out = DirectSum({out gamma_local});
end

end


