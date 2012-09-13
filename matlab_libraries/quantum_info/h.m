function[out] = h(x,option);
if (option=='bosons')
	if(x==1)
		out=0;
	else
		out= ((x+1)/2)*log2((x+1)/2) - ((x-1)/2)*log2((x-1)/2);
	end
else 
	if(x==0)
		out=0;
	else
		out=-x*log2(x);	
	end
end
