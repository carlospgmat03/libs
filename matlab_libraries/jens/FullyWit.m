function[c Z] = FullyWit(gamma,n,constraints)
% function[c Z]=FullyWit(gamma,n,constraints)
%
% Written by Philipp Hyllus and Jens Eisert in September 2005
% Missing function 'sigma' added in March 2006
%
% This function tests whether the multipartite covariance matrix 
% gamma is fully separable or not. If it is not, then the output 
% matrix Z is an entanglement witness that can be used to 
% prove this fact. 
%
% Input:    'gamma' is the covariance matrix of a state
%           that is held by the parties. 
%
%           'n' is a vector with the numbers of modes that the parties hold.
%           If A holds 3 modes, B 2 modes, and C 4 modes, then n=[3 2 4].
%
%           'Constraints' gives extra experimental constraints
%           of the form Tr[Z*A]=0 on the witness Z. 
%           If no extra constraints are wanted, then constraints=0.
%           Otherwise, the constraint matrices A1,A2,...,Ak
%           have to be written into constraint as [A1 A2 ... Ak].
%           Note: there might be errors if the matrices Ai are defined as
%           sparse!
%
% Output:   c=tr[gamma Z]-1. If c>=0 then gamma
%           is separable and entangled otherwise.
%           The entanglement witness Z is also returned.

% Initialization
N=sum(n);  % Total number of modes
P=size(n,2); % Number of parties

Z1=sdpvar(2*N,2*N,'hermitian','complex');
Z2=sdpvar(2*N,2*N,'hermitian','complex');

% Constraints
F=set(Z1>0)+set(Z2>0);
F=F+set(trace(Z2*i*sigma(N))==-1);

% Equality of block diagonal real parts
j=1;
for k=1:P,
    F=F+set(real(Z1(j:j+2*n(k)-1,j:j+2*n(k)-1))==real(Z2(j:j+2*n(k)-1,j:j+2*n(k)-1)));
    j=j+2*n(k);
end

% Additional constraints
if size(constraints,2)>1
    for k=1:(size(constraints,2)/(2*N)),
        F=F+set(trace(real(Z1))*constraints(:,2*N*(k-1)+1:2*N*k)==0);
    end
end

% Show the constraints
F
pause

%Go!
solvesdp(F,trace(Z1*gamma)+trace(Z2*i*sigma(N)))

% output
c=double( (trace(Z1*gamma)+trace(Z2*i*sigma(N))))
Z=real(double(Z1))

%==================================================================
function[out]=sigma(N)
% gives the sigma matrix for N modes
out=sparse(2*N,2*N);
for k=1:N,
    out(2*k-1,2*k)=1;
    out(2*k,2*k-1)=-1;
end