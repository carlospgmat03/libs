function[out] = DistanceFromSeparableStates(gamma, modesAB, type);
% Function to calculate entanglement. The calculation is based on the idea 
% by Jens 
% http://www.iop.org/EJ/article/1367-2630/8/4/051/njp6_4_051.html
% We kinda propose a SDP problem, and then solve it using Matlab routines. 
% Determin the number of modes. 

N = sum (modesAB);
if type ==1 
	Gamma1=sdpvar(N,N,'symmetric');
	Gamma2=sdpvar(N,N,'symmetric');
	xe=sdpvar(1,1);
	F = set([]);
	F=F + set(gamma - [Gamma1 zeros(N); zeros(N) Gamma2]   > 0);
	F=F + set(Gamma1 + i*qosigma(N/2) > xe*eye(N));
	F=F + set(Gamma2 + i*qosigma(N/2) > xe*eye(N));
	solvesdp(F, -xe , sdpsettings('verbose',0, 'debug',1)   );
	out = double(xe); 
elseif type == 2
	Gamma1=sdpvar(N,N,'symmetric');
	Gamma2=sdpvar(N,N,'symmetric');
	xe=sdpvar(1,1);
	% the constrains. see (31) of the paper. 
	F = set([]);
	F=F + set(gamma - [Gamma1 zeros(N); zeros(N) Gamma2]   > 0);
	F=F + set(Gamma1 + (1+xe)*i*qosigma(N/2) > 0);
	F=F + set(Gamma2 + (1+xe)*i*qosigma(N/2) > 0);
	solvesdp(F, -xe , sdpsettings('verbose',0, 'debug',1)   );
	out = double(xe); 
end 
