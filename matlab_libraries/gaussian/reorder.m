function[out] = reorder(Gamma, modes);
totalmodes=length(Gamma)/2;
if (size(modes) ~= totalmodes)
  error('Modes especification do not compile. Not enough modes!');
end
if ( norm(sort(modes)-(1:totalmodes)) >0)
  error('Modes especification do not compile. Not unique modes. ');
end
y=[2*modes-1 ; 2*modes];
y=y(:)';
out=Gamma(y,y);
return
