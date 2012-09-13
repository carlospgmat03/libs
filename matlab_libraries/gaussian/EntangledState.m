function[out]= EntangledState(x)

alpha = x/(1-x);
s = (exp(alpha) - exp(-alpha) )/2 ;
c = (exp(alpha) + exp(-alpha) )/2 ;
out = [c 0 s 0; 0 c 0 -s;  s 0 c 0;   0 -s 0 c ];
%out = c^2 - s^2;
end

