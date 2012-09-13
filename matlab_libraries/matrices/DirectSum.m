function[out]= DirectSum(A)
%function[out]= DirectSum(A)
%Direct sum of several matrices, Say 
%A = {g g g2} 
%DirectSum(A) 
%produces the direct product of g, g, and g2

out = [];
for k=1:length(A)
	[OldRows, OldCols] =  size(out);
	[NewRows, NewCols] = size(A{k});
	out = [out zeros(OldRows, NewCols); zeros(NewRows, OldCols) A{k}];
end
