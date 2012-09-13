function[out] = IsEntangled(gamma, modesAB, epsilon);

na = modesAB(1);
nb = modesAB(2);
xe=sdpvar(1,1,'symmetric');
Gamma1=sdpvar(2*na,2*na,'symmetric');
Gamma2=sdpvar(2*nb,2*nb,'symmetric');
F = set([]);
F=F + set(gamma - [Gamma1 zeros(2*na, 2*nb); zeros(2*nb, 2*na) Gamma2]   > 0);
F=F + set(Gamma1 + i*xe*qosigma(na) > 0);
F=F + set(Gamma2 + i*xe*qosigma(nb) > 0);
sol = solvesdp(F, -xe , sdpsettings( 'verbose',0)   );


out = (double(xe) < 1-epsilon);
end



