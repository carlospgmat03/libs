% I want to calculate the stat that is farthest both from the 
% separable states and from the distilable states. We want to state it as a SDP

clear; 
path(path, 'cfp/directories');
all_cfp_dirs

% First the two distances that we want to use. 
% Here we assume that each of the two parties has two modes. 
modesAB=[2 2] ;
% modesAB=[1 1] ;

fprintf('We want to find a state with some distance from distilable states. \n');
fprintf('The state must have more than ammount of entanglement. \n');
AmmountOfEntanglement=10.1;
N = sum (modesAB);
% The variables we want to diagonaliza over:
xe = AmmountOfEntanglement;
gamma =sdpvar(2*N,2*N,'symmetric');
Gamma1=sdpvar(  N,  N,'symmetric');
Gamma2=sdpvar(  N,  N,'symmetric');

F = set([]);
F=F + set(gamma - [Gamma1 zeros(N); zeros(N) Gamma2]   > 0);
F=F + set(Gamma1 + (1+xe)*i*qosigma(N/2) > 0);
F=F + set(Gamma2 + (1+xe)*i*qosigma(N/2) > 0);

solvesdp(F, -xe , sdpsettings('verbose',0, 'debug',1)   );
gamma


return;

fprintf('We want to find a state with more than some given entanglement. \n');
fprintf('The state must have more than ammount of entanglement. \n');
AmmountOfEntanglement=10.1;
N = sum (modesAB);
% The variables we want to diagonaliza over:
xe = AmmountOfEntanglement;
gamma =sdpvar(2*N,2*N,'symmetric');
Gamma1=sdpvar(  N,  N,'symmetric');
Gamma2=sdpvar(  N,  N,'symmetric');

F = set([]);
F=F + set(gamma - [Gamma1 zeros(N); zeros(N) Gamma2]   > 0);
F=F + set(Gamma1 + (1+xe)*i*qosigma(N/2) > 0);
F=F + set(Gamma2 + (1+xe)*i*qosigma(N/2) > 0);

solvesdp(F, -xe , sdpsettings('verbose',0, 'debug',1)   );
gamma


return;

N = sum (modesAB);
gamma=sdpvar(2*N,2*N,'symmetric');
F = set([]);
F=F + set(gamma + i*qosigma(N) > 0);

% Distance 1. the distance from the separable states:
% 
% 
Gamma_separability_1=sdpvar(N,N,'symmetric');
Gamma_separability_2=sdpvar(N,N,'symmetric');
% distance_separability=sdpvar(1,1);
distance=sdpvar(1,1);
% the conditions to find separability
F=F + set(gamma - [Gamma_separability_1 zeros(N); zeros(N) Gamma_separability_2]   > distance*eye(2*N));
% F=F + set(Gamma_separability_1 + (1-distance)*i*qosigma(N/2) > 0);
% F=F + set(Gamma_separability_2 + (1-distance)*i*qosigma(N/2) > 0);
F=F + set(Gamma_separability_1 + i*qosigma(N/2) > 0 );
F=F + set(Gamma_separability_2 + i*qosigma(N/2) > 0 );
% quremos que distance sea como 

% Distance 2, the distance from distilable states:
m_distance=sdpvar(1,1);
F=F + set(gamma - i* PPT(qosigma(N), modesAB ) > m_distance * eye(2*N)  );
F=F + set(gamma - i* PPT(qosigma(N), modesAB ) > m_distance * eye(2*N)  );

t_distance=sdpvar(1,1);
F=F + set([-distance 0;  0 -m_distance ]  > t_distance * eye(2)  );
% Tambien falta poner qeu la matriz sea fisica

% solvesdp(F, -m_distance   );d1=double(m_distance);
% F


%solvesdp(F, -distance   ); d1=double(distance)
solvesdp(F,  - t_distance , sdpsettings('verbose',0, 'debug',1)   )


%overall_distance=sdpvar(1,1);
%F=F + set( [m_distance > zeros(N); zeros(N) Gamma_separability_2]  > overall_distance*eye(2));

% creo que esto viendo el problema de min max

g=double(gamma)
 

d_s = DistanceFromSeparableStates(g, modesAB)
fprintf('distancia de separables: %e \n', d_s);
% fprintf('distancia de distilable: %e \n', DistanceFromDistilableStates(g, modesAB));

