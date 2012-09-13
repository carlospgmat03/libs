function[c Z] = MultiWit(gamma,n,constraints)
% function[c Z] = MultiWit(gamma,n,constraints)
%
% Written by Philipp Hyllus and Jens Eisert in September 2005.
% Performance substantially improved with the help of 
%   Johan Loefberg in March 2006.
%
% This function tests whether the multipartite covariance matrix 
% gamma is biseparable or not. If it is not, then gamma is
% truly multipartite entangled and the output matrix Z is an 
% entanglement witness that can be used to prove this fact. 
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
% Output:   c=tr[gamma*Z]-1. If c>=0 then gamma
%           is biseparable and multiparty entangled otherwise.
%           The entanglement witness Z is also returned.

%Initialization
N=sum(n);       % Total number of modes
P=size(n,2);    % Number of parties
K=2^(P-1)-1;    % Total number of distinct bipartitions

Z1=sdpvar(2*N*(K+1),2*N*(K+1),'hermitian','complex'); 
% The matrix holding the matrices.
%
% Note: This parametrization leads obviously to a waste of space, 
% we are currently working on a better one.

ZS=sdpvar(1,2+K); % this vector holds the scalars needed.

% The following loops sets the off-block elements to zero.
for j=1:K,   
 Z1(2*N*(j-1)+1:2*N*j,2*N*j+1:2*N*(K+1))=0;
 Z1(2*N*j+1:2*N*(K+1),2*N*(j-1)+1:2*N*j)=0;   
end                                                 
F = unblkdiag(set(Z1>0));

F=F+set(ZS>0);
F=F+set(ZS(1)-ZS(2)==1); 

k=1;
for j=1:K,
    F=F+set(trace(i*sigma(N)*Z1(2*N*j+1:2*N*(j+1),2*N*j+1:2*N*(j+1)))+ZS(1)-ZS(2)+ZS(2+j)==0);

    % Equality of block diagonal real parts of all bipartite partitions
    F=F+set(real(CutParties(Z1(1:2*N,1:2*N),k))==real(CutParties(Z1(2*N*j+1:2*N*(j+1),2*N*j+1:2*N*(j+1)),k)));
    F=F+set(real(CutParties(Z1(1:2*N,1:2*N),invertk(k,P)))==real(CutParties(Z1(2*N*j+1:2*N*(j+1),2*N*j+1:2*N*(j+1)),invertk(k,P))));
    k=inck(k,P);
end

% Additional constraints
if size(constraints,2)>1
    for k=1:(size(constraints,2)/(2*N)),
        F=F+set(trace(real(Z1(1:2*N,1:2*N))*constraints(:,2*N*(k-1)+1:2*N*k))==0);
    end
end

F
pause

% Go!
solvesdp(F,trace(Z1(1:2*N,1:2*N)*gamma)-ZS(1)+ZS(2))

% Output
c=double(trace(Z1(1:2*N,1:2*N)*gamma)-ZS(1)+ZS(2));
Z=double(real(Z1(1:2*N,1:2*N)))


%-------------- further functions -----------------------------------------
function[out]=sigma(N)
% gives the sigma matrix for N modes
out=sparse(2*N,2*N);
for k=1:N,
    out(2*k-1,2*k)=1;
    out(2*k,2*k-1)=-1;
end

function[out]=CutParties(gamma,n);
% function[out]=CutParties(gamma,n)
%
% Inputs: 'gamma' is a covariance matrix
%         'n'     is a vector of integers
% Output: This function returns
%         a smaller covariance matrix, where
%         the rows and columns of the parties
%         in 'n' are cut.
%
% See also: CutMat

for k=1:size(n,2),
    v(2*k-1)=2*n(k)-1;
    v(2*k)=2*n(k);
end
out=CutMat(gamma,v);

%---------------------------------------------------------
function[out]=CutMat(M,n);
% Inputs: M is a matrix
%         n is a vector of integers
%
% Description: This function returns
% a smaller hermitian matrix, where
% the rows and columns in n are cut.
% 
% Example: M=[1 2 3; 2 1 4; 3 4 1]
%          n=[1 2]
% CutMat(M,n)=[1]

help=M;
N=size(n,2);
out=CutCol(CutRow(help,n(1)),n(1));
if N>1
    out=CutMat(out,n(2:N)-1);
end

%---------------------------------------------------------
function[out]=CutCol(M,n);
% function[out]=CutCol(M,n);
%
% Inputs: M is a matrix
%         n is an integer
%
% Description: This function returns
% a matrix where the n'th column is
% taken out.

r=size(M,1);
c=size(M,2);
help=M;
out=[help(1:r,1:(n-1)) help(1:r,(n+1):c)];

%---------------------------------------------------------
function[out]=CutRow(M,n);
% function[out]=CutRow(M,n);
%
% Inputs: M is a matrix
%         n is an integer
%
% Description: This function returns
% a matrix where the n'th row is
% taken out.

r=size(M,1);
c=size(M,2);
help=M;
out=[help(1:(n-1),1:c); help((n+1):r,1:c)];

%---------------------------------------------------------
function[out]=inck(k,n)
% function[out]=inck(k,n)
%
% this function gives the next partition
% after the one carried by the vector k
% if k is the maximal vector already,
% and L<n, then out(1)=[1 2 3 ... L+1] where 
% L is the length of k.
% 
% Examples: inck([1 2],5)=[1 3]
%           inck([2 5 6],6)=[3 4 5]
%           inck([4 5],5)=[1 2 3]

L=size(k,2);
if (k(1)>(n-L)) %then k is the last vector
    if (L<n)
        out=incvec(L+1,1);
    else
        out=[n k(2:L)]; 
    end
else
    j=L;
    while k(j)>(n-(L-j+1))
        j=j-1;
    end
    out=[k(1:(j-1)) incvec(L-j+1,k(j)+1)];
end

%---------------------------------------------------------
function[out]=incvec(l,s)
% function[out]=incvec(l,s)
%
% returns a vector of length 'l'
% which entries start with 's' and increase by 1.
out=zeros(1,l);
for j=1:l,
    out(j)=s+j-1;
end

%---------------------------------------------------------
function[out]=invertk(k,n)
% function[out]=invertk(k,n)
% Input: 'n' is the total number of parties
%        'k' is a vector holding a subset of parties
%            should be increasingly ordered
% Output: the vector where all parties are collected which are not in 'n'

% first the vector holds all parties
out=zeros(1,n);
for j=1:n,
    out(j)=j;
end

% and then the parties in 'k' are cut out
for j=1:size(k,2),
    out=[out(1:(k(j)-1)) out((k(j)+1):(n-j+1))];
    k=k-1;
end 