function[out]= OperatorOnVariables(Operator, ModesToPutAt, N)
% 
rowscols = zeros(1, 2*length(ModesToPutAt));
rowscols(1:2:length(rowscols)) = 2*ModesToPutAt - 1;
rowscols(2:2:length(rowscols)) = 2*ModesToPutAt ;

out = eye(2*N);
out(rowscols, rowscols) = Operator;

end
