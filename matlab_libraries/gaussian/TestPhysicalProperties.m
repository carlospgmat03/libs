% In this script we test numerically that the states have the 
% desired Phyiscal properties. All parameters should be between 0 and 1. 



% Here we test that the operator fullfiles 
% \gamma - \rmi \sigma \ge 0

% The thermal states, which are mutiples of the identity. 
% PhysicalityCorrelationMatrices(EntangledState(1.3))
clear ;
path(path, 'cfp/directories');
all_cfp_dirs
close all;
fplot ('PhysicalityCorrelationMatrices(ThermalState(z,4))',[0,.95],'b');
hold on;
%fplot ('PhysicalityCorrelationMatrices(EntangledState(z))',[0,.95],'LineWidth',2.);
fplot ('PhysicalityCorrelationMatrices(EntangledState(z))',[0,.95],'r');
title('Minimum eigenvalue of matrix \gamma - i \sigma');
legend('Thermal states', 'Entangled States');

% Now an entanglement test for the states. 

figure(2);
hold on;
fplot ('DistanceFromDistilableStates(ThermalState(z, 4), [1 1])',[0.01,.9],'b');
fplot ('DistanceFromDistilableStates(EntangledState(z), [1 1])',[0.01,.95],'r');
title('Distance From Distilable States');
%axis([0 1 -.5 1]);
legend('Thermal states', 'Entangled States');


figure(3);
hold on;
fplot ('DistanceFromSeparableStates(ThermalState(z, 4), [1 1])',[0.01,.9],'b');
fplot ('DistanceFromSeparableStates(EntangledState(z), [1 1])',[0.01,.95],'r');
title('Distance From Separable States:');
axis([0 1 -.5 1]);
legend('Thermal states', 'Entangled States');


