function Gamma = CorrelationMatrix(a,StringSpecification)

mr1=[a(9) a(10) ; a(11) a(12) ];
mr2=[a(13) a(14) ; a(15) a(16) ];
dr1=diag([a(5) a(6)]);
dr2=diag([a(7) a(8)]);
	
Gamma = [a(1)*eye(2) 	zeros(2)  	dr1 		mr1; ...
	zeros(2) 	a(2)*eye(2)  	mr2 		dr2; ...
	dr1' 		mr2' 		a(3)*eye(2) 	zeros(2); ...
	mr1' 		dr2' 		zeros(2) 	a(4)*eye(2) ];
if (nargin>1 & StringSpecification=='Least possible noise')
	Gamma = - min(real(eig(Gamma + i*(qosigma(4) )) ))* eye(8) + Gamma;
end


