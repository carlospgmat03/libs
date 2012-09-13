function [A, S] = standard_form_2_modes(Gamma, debug);
% J. Eisert, C. Pineda Nov. 6, 2009
if (nargin<2)
	debug = false;
end 


if debug; fprintf('To 2 mode standard form\n'); end

A=Gamma;
S=eye(4);
z2=zeros(2,2);

V1=orthogonal_2_2_diagonalizator(A(1:2,1:2));
% [V1, dummy]=eig(A(1:2,1:2))
% if (V1(1,1)*V1(2,2)<0); V1(1,:)=-V1(1,:); end
% [V2, dummy]=eig(A(3:4,3:4));
% if (V2(1,1)*V2(2,2)<0); V2(1,:)=-V2(1,:); end
V2=orthogonal_2_2_diagonalizator(A(3:4,3:4));
S_step=[V1 z2; z2 V2];
S= S*S_step;
A=S_step'*A*S_step;

d1=(A(2,2)/A(1,1))^(1/4);
d2=(A(4,4)/A(3,3))^(1/4);
S_step=[diag([d1, 1/d1]) z2; z2 diag([d2 1/d2])]; S= S*S_step;
A=S_step'*A*S_step;
if debug;
	fprintf('Norm of symplectic matrix: %f\n',...
	norm(S'*qosigma(2)*S-qosigma(2))); 
end
A1=A;
B=A(1:2,3:4);
[U, dummy, V]= svd(B);
if (U(1,1)*U(2,2)<0); U(:,1)=-U(:,1); end
if (V(1,1)*V(2,2)<0); V(:,1)=-V(:,1); end
S_step  = [U z2; z2 V]; S= S*S_step;
A=S_step'*A*S_step;
distance_to_symplectic(S_step);


function O=orthogonal_2_2_diagonalizator(b)
	% devuelve O t.q. O'*b*O es diagonal
	x = 2*b(1,2)/(b(1,1)-b(2,2));
	if (abs(x) < 1)
		theta = 0.5*atan(x);
	else
		theta= 0.5*acot(0.5*(b(1,1)-b(2,2))/b(1,2));
	end
	ct=cos(theta);
	st=sin(theta);
	O=[ct -st; st ct];

function d = distance_to_symplectic(S);
	n = size(S,2)/2;
	%mysigma = [ zeros(n,n) eye(n); -eye(n) zeros(n,n)];
	%d=norm(S'*mysigma*S-mysigma);
	d=norm(S'*qosigma(n)*S-qosigma(n));
 
