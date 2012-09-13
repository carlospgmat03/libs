function[out] = Translate(Data, KindOfTranslation)
if (KindOfTranslation == 'TransmitanceSqueezingFromSqueezingNoise')
  d = Data(1); % Pure Squeezing
  a = Data(2); % Pure Noise
  Squeezing = (d-d^2-d*a)/(1-d+d*a);
  Transmitance = (1-2*d+d^2-a+2*d*a-d^2*a-d*a^2)/(1-2*d+d^2+2*d*a);
  out = [Transmitance Squeezing];
end
