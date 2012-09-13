function [K,A,N] = iFactor(S); 
% M. Benzi, N. Razouk / Applied Mathematics Letters 20 (2007) 260–265 265 
% 
% This function computes the Iwasawa decomposition of 
% a real symplectic matrix of order 2n. 
% 
% Input: a real symplectic matrix [S_11 S_12; S_21 S_22] 
% 
% Output: K = 2n-by-2n orthogonal symplectic matrix 
% A = 2n-by-2n positive diagonal symplectic matrix 
% N = 2n-by-2n "triangular" symplectic matrix 
% 
% s.t. 
% 
% S = K*A*N 
% 
n_2 = size(S); 
n = n_2/2; 
% Compute thin QR factorization of S1 = [S_11; S_21]. 
S_11 = S(1:n,1:n); 
S_21 = S(n+1:n_2,1:n); 
S1 = [S_11; S_21]; 
[Q,R] = qr(S1,0); 
Q_11 = Q(1:n,1:n); 
Q_21 = Q(n+1:n_2,1:n); 
% Compute U and D from given R where U is unit upper triangular 
% and H is a diagonal matrix such that R = H*U and 
% R’*R = U’*H^2*U. 
H = diag(diag(R)); 
U = H\R; 
% Compute blocks for the factors K, A, N. 
h = sign(diag(H)); 
SQRT_D = diag(h.*diag(H)); 
SQRT_D_inv = diag(1./diag(SQRT_D)); 
h = h'; 
h = h(ones(1,n),:);
K_11 = Q_11.*h; 
K_12 = -Q_21.*h; 
% Form the Iwasawa factors K, A, N. 
A = [SQRT_D zeros(n) ; zeros(n) SQRT_D_inv]; 
K = [K_11 K_12; -K_12 K_11]; 
S1 = S(1:n_2,n+1:n_2); 
N1 = A\(K'*S1); 
N = [U N1(1:n,1:n); zeros(n) N1(n+1:n_2,1:n)]; 

