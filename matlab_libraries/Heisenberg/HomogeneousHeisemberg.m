function[ham] = HomogeneousHeisemberg(J,N)

%HomogeneousHeisemberg(J, N): Creation of the homogeneous open Heisenberg
%chain.

ham=zeros(2.^N,2.^N);

for i=0:N-2
    
    ham=J.*Heisenbergterm(i,N)+ham;
    
end

