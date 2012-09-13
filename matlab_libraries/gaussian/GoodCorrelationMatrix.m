function[out]= GoodCorrelationMatrix(gamma, epsilon)
N=size(gamma,1);
out = min(eig( gamma + i*qosigma(N/2) ) ) ;
end



