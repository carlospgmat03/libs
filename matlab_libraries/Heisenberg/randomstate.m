function[state] = randomstate(N)

%randomstate(N): Random quantum state of dimension N.

    state=rand([N,1])+1j*rand([N,1]);
    
    normalization=norm(state);
    
    state=normalization.^(-1).*state;