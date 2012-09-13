function[out] = RandomCorrelationMatrix(Modes);
Gamma=rand(2*Modes)-0.5; 
Gamma=Gamma+Gamma';
out=Gamma+eye(2*Modes)*(-min(eig(Gamma-i*qosigma(Modes)))+rand());

