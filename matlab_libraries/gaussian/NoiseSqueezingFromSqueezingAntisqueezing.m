function[out] = NoiseSqueezingFromSqueezingAntisqueezing(Squeezing, AntiSqueezing);


S=Squeezing;
A=AntiSqueezing;

x=power(10,S/10)-power(10,A/10);
d=(x+sqrt(x^2+4))/2;
alpha=power(10,S/10)-d;
out=[d alpha];

