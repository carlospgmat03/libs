function Gam=BoundEntPPT

c=rand(1);
e=rand(1);
b=rand(1);
f=rand(1);
aa=rand(1);
a=cos(aa)^2*c*e;

Gam=[(b+e*f/(-a+c*e))/a 0 0 0 -f/(a-c*e) 0 0 0; 0 a/b 0 0 0 0 0 -1/b; 0 0 (b+e*f/(-a+c*e))/a 0 0 0 f/(a-c*e) 0; 0 0 0 a/b 0 -1/b 0 0 ; -f/(a-c*e) 0 0 0 c*f/(-a+c*e) 0 0 0 ; 0 0 0 -1/b 0 (b*c+a*b*e+f)/(a*b*f) 0 0; 0 0 f/(a-c*e) 0 0 0 c*f/(-a+c*e) 0;  0 -1/b 0 0 0 0 0 (b*c+a*b*e+f)/(a*b*f)];
end
