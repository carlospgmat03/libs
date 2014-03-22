function [mat] = paulimatrices(n)

%paulimatrices(i) indexed like mathematica, but without paulimatrices(0).

switch n
    
    case 3
        mat=[1,0;0,-1];
        
    case 1
        
        mat=[0,1;1,0];
        
    case 2
        
        mat=[0,-1i;1i,0];
        
end