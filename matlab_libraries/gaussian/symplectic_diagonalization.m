function [S, EigenValues] = symplectic_diagonalization(Gamma, debug);
% J. Eisert, C. Pineda Nov. 6, 2009
% I think one must follow section III of J. Math. Phys., Vol 40, No.7, July 1999
% by Simon, Chaturvedi, and Srinivasan
% So its the same principle but the "simplectic identity" induced 
% by the order (x1 p1 x2 p2 ... xn pn)
% Dependencies: qosigma, that provides the skew-symmetric matrix
% corresponding to the mentioned order



if (nargin<2); debug = false; end 

n = size(Gamma,2)/2;
cumulative_error=0;
beta=qosigma(n);

% Build the antisymmetric "M" matrix that does the job
V_m12=inv(sqrtm(Gamma));
M=(V_m12*beta*V_m12); 
%Thus this eigenvalues come in pairs with 
% i lambda, -i lambda
[R, D]=eig(M) ;
D=diag(D);
cumulative_error=cumulative_error + norm(real(D));
cumulative_error=cumulative_error + norm(D(1:2:end)+D(2:2:end));
% so it is a matter of reorganization to set up things in
% the right order.


% This block diagonal matrix reorders the eigenvalues 
% in such a way as to have the correct order.
joto=[];
for j1=1:n; joto=two_block_diagonal(joto,[1, i; i 1]); end
S=V_m12*R*joto*(-1)^(-1/4);
cumulative_error=cumulative_error + norm(imag(S));
S=real(S);

% Get the actual values that correspond to the 
% eigenvalues. 
C=S'*beta*S;
cumulative_error=cumulative_error + norm(imag(C));
C=real(C);
D=diag(kron(1./sqrt(diag(C(2:2:2*n,1:2:2*n))),[1; -1]));
S=S*D;



% Verify that the state is indeed diagonalized by the 
% matrix and that it is indeed symplectic
Test=S'*Gamma*S;
EigenValues=diag(Test);
cumulative_error=cumulative_error + norm(Test-diag(EigenValues));
cumulative_error=cumulative_error + distance_to_symplectic(S);

% Test that the eigenvalues are the same as with the good old method
eigenvalues = real(eig(i*beta*Gamma));
eigenvalues = kron(sort(eigenvalues(eigenvalues >0)),[1 ; 1]);
% eigenvalues = kron(eigenvalues(),[1 ; 1]);
cumulative_error = cumulative_error + norm(sort(diag(Test))-sort(eigenvalues));

% S
% EigenValues

if (cumulative_error > 10e-13)
	error(' symplectic_diagonalization:BadDiag', ...
		'There is a problem. Some error')
end



function C=two_block_diagonal(A,B)
	C=[A zeros(size(A,1), size(B,2) ); zeros(size(B,1),size(A,2)) B];
function d = distance_to_symplectic(S);
	n = size(S,2)/2;
	mysigma = [ zeros(n,n) eye(n); -eye(n) zeros(n,n)];
	d=norm(S'*mysigma*S-mysigma);
	d=norm(S'*qosigma(n)*S-qosigma(n));

