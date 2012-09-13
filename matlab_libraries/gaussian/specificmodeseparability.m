function[out] = specificmodeseparability(Gamma, modes);
% Check the value of entanglement for a given bi-separability 
% of the modes indicated in the array modes
totalmodes=length(Gamma)/2;
if (strcmp(class(modes),'cell'))
 % First, join all the modes in an array. Check then than it is
 % the same as an unorder set range totalmodes
 neworder = concatenate(modes);
 reorderedGamma = reorder(Gamma, neworder);

if exist('RandStream')
	previous_random_stream = RandStream.getDefaultStream();
end
F = set([]);
number_of_parties=length(modes);
Gsep=cell(number_of_parties,1);
for j=1:number_of_parties
  modes_j = length(modes{j});
  Gsep(j)={sdpvar(2*modes_j,2*modes_j,'symmetric')};
end
GammaSeparable =  DirectSum(Gsep);
F=F + set(reorderedGamma - GammaSeparable   > 0);
xe=sdpvar(1,1,'symmetric');
F=F + set(GammaSeparable + i*xe*qosigma(totalmodes) > 0);
sol = solvesdp(F, -xe , sdpsettings( 'verbose',0));
out = double(xe);
if exist('RandStream')
	RandStream.setDefaultStream(previous_random_stream);
end


 return
end 
modesother=setdiff([1:totalmodes], modes);
y=[(2*modes-1) (2*modesother-1) ; 2*modes 2*modesother];
y=y(:)';
Gamma2=Gamma(y,y);
out = how_entangled(Gamma2, [length(modes) length(modesother)]);
return



