function[out] = concatenate(data)

if (strcmp(class(data),'cell'))
 nn=length(data);
 out =[];
 for j=1:nn
  out=[out data{j}];
 end
 return

end
