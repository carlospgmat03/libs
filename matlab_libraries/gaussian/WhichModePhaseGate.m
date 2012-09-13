function[out]= WhichModePhaseGate(phi, k, N)
% 
out = DirectSum({eye(2*(k-1)) SingleModePhaseGate(phi) eye(2*(N-k))});
end
