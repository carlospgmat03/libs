function[out] = PPTnessMeasureCorrelationMatrices(gamma, modesAB, type);
% Function to calculate entanglement. The calculation is based on the idea 
% by Jens 
% http://www.iop.org/EJ/article/1367-2630/8/4/051/njp6_4_051.html
% We kinda propose a SDP problem, and then solve it using Matlab routines. 
% Determin the number of modes. 
N = sum (modesAB);
% The variables we want to diagonaliza over:
%Gamma1=sdpvar(N,N,'symmetric');
%Gamma2=sdpvar(N,N,'symmetric');

if type ==1 
	xe=sdpvar(1,1);
	% the constrains. see (31) of the paper. 
	F = set([]);
	F=F + set(gamma + i* PPT(qosigma(N), modesAB ) > xe * eye(2*N)  );
	solvesdp(F, -xe , sdpsettings('verbose',0, 'debug',1)   );
	% out = max([0 +double(xe)]); 
	out = double(-xe); 
elseif type ==2 
	xe=sdpvar(1,1);
	% the constrains. see (31) of the paper. 
	F = set([]);
	F=F + set(gamma + i*(1+xe)* PPT(qosigma(N), modesAB ) > 0  );
	solvesdp(F, -xe , sdpsettings('verbose',0, 'debug',1)   );
	% out = max([0 +double(xe)]); 
	out = double(-xe); 
end 
