function[out] = ApplyPhase(Gamma, phi, n);

NumberOfModes = size(Gamma,2)/2;
Gate = WhichModePhaseGate(phi, n, NumberOfModes);
out= Gate*Gamma*Gate';
end


