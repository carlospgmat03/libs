function[out] = SqueezingAntisqueezingFromdalpha(d, alpha);


a1=10*log10(d+alpha);
% 
% alpha
% d
% (1./d)+alpha
a2=10*log10((1./d)+alpha);

Squeezing=min(a1,a2);
AntiSqueezing=max(a1,a2);

out=[Squeezing AntiSqueezing];


