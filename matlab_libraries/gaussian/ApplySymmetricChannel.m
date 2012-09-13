function[out] = ApplySymmetricChannel(Gamma, phases);

ModesSys=size(phases,2);

out=Gamma;
for j1=1:ModesSys
	out=ApplyBeamSplitter(out, [j1 ModesSys+j1], 0, pi/4, phases(j1));
end

end


