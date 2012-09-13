function[out]= ThermalState(alpha, N)
% alpha should be between 0 and 1. If not, return error. 
% 

beta = 1./(1.*(1-alpha));
out = beta.* eye(N);
end

