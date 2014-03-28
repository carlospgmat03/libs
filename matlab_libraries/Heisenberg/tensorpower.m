function [newstate] = tensorpower(state,n)

%tensorpower(state,n)
%where n is the exponent.

newstate=state;

for i=0:n-2
    
    newstate=kron(newstate,state);

end
