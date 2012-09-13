function [eigenvalues, S, errorvalue] = symplectic_2_mode_symmetric_diagonalization(Gamma, debug);
% J. Eisert, C. Pineda Nov. 6, 2009
% 2 mode following from JPB, 37 (2004) L21
% Symplectic invariants, entropic measures and correlations of Gaussian states
% Alession Serafini, and others. 
if (nargin<2)
	debug = false;
end 

% Step zero, reduce to standard form
[A, S] = standard_form_2_modes(Gamma, debug);

% First a rotation to eliminate the (1,3) element
phi=0.5*atan2(2*A(3,1),A(1,1)-A(3,3));
S_step=two_mode_rotation(phi);
A=S_step'*A*S_step;
S=S*S_step;

% Then a two squeezing to equate elements (1,1) and (3,3)
r=(A(3,3)/A(1,1))^(1/4);
S_step=Squeezing_two_mode(r);
A=S_step'*A*S_step;
S=S*S_step;

% Finally a rotation to eliminate the (2,4) element
phi=0.5*atan2(2*A(4,2),A(2,2)-A(4,4));
S_step=two_mode_rotation(phi);
A=S_step'*A*S_step;
S=S*S_step;

% We rotate with local squeezing to get a purely thermal state
S_step = Squeezing([(A(2,2)/A(1,1))^(1/4) (A(4,4)/A(3,3))^(1/4)]);
A=S_step'*A*S_step;
S=S*S_step;

eigenvalues = diag(A);
errorvalue = distance_to_symplectic(S)+norm(A-diag(eigenvalues)) +norm(S'*Gamma*S-A);

function R=two_mode_rotation(phi)
	c=cos(phi); s=sin(phi);
	R=[c 0 -s 0; 0 c 0 -s; s 0 c 0; 0 s 0 c];

function Stm=Squeezing(r)
	Stm=[ ];
	for j1=1:size(r,2)
		Stm=[Stm [r(j1) 1/r(j1)]];
	end
 	Stm=diag(Stm);


function Stm=Squeezing_two_mode(r)
 	Stm=diag([ r 1/r 1/r r ]);

function d = distance_to_symplectic(S);
	n = size(S,2)/2;
	%mysigma = [ zeros(n,n) eye(n); -eye(n) zeros(n,n)];
	%d=norm(S'*mysigma*S-mysigma);
	d=norm(S'*qosigma(n)*S-qosigma(n));


% function C=block_diagonal(A)
% 	for int 
% 	C=[A zeros(size(A,1), size(B,2) ); zeros(size(B,1),size(A,2)) B];
function C=two_block_diagonal(A,B)
	C=[A zeros(size(A,1), size(B,2) ); zeros(size(B,1),size(A,2)) B];


