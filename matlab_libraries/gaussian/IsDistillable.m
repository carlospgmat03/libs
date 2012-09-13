function[out] = IsDistillable(gamma, modesAB, epsilon);

N = sum (modesAB);
%eig( gamma + i* PPT(qosigma(N), modesAB ))  
out = (min (real(eig( gamma + i* PPT(qosigma(N), modesAB ))  )) < -epsilon);
%out = min(eig( gamma + i* PPT(qosigma(N), modesAB )))  ;
%out = (  min(eig( gamma + i* PPT(qosigma(N), modesAB )))  < -epsilon );
end
 
