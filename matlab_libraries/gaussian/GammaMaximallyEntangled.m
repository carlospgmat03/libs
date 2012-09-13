function Gamma=GammaMaximallyEntangled(n,r);
% Maximally entangled acording to PRA 66, 032316 Eq. (4)
% r=4.;
% n=2;

Lambda  = diag(kron(ones(n,1)', [1 -1]));
Ar = cosh(r)* eye(2*n);
Cr = sinh(r*Lambda);
Gamma = [ Ar Cr; Cr Ar];

