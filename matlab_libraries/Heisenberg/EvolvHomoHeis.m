function[finalstate] = EvolvHomoHeis(J,state,t)

%EvolvHomoHeis(J, state, t): Evolution of the homogeneous Heisemberg open
%chain

    finalstate=state;
    
    qubits=log2(length(state));
    
    H=HomogeneousHeisemberg(J,qubits);
    
    U=expm(-1j.*H.*t);
    
    finalstate=U*finalstate;