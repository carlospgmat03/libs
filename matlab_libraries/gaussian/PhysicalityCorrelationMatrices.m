function[out]=PhysicalityCorrelationMatrices(gamma)


NumberOfModes = size(gamma,1)/2;
mine = min (eig(gamma + i*qosigma(NumberOfModes))) ; 
out = (mine > 0);
end


