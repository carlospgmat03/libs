function [term] = Heisenbergterm(i,N)

%Heisenbergterm(i,N): Calula \vec sigma_i \cdot \vec sigma_{i+1} donde N es
%el numero de qubits.

term=0;

    if N<2+i
        
       'indice de termino muy grande, procure que sea i<=N-2'
       
    end

    if N>2+i
        
        C=tensorpower(eye(2,2),N-2-i);
        
        for j=1:3
            
           term=kron(tensorpower(paulimatrices(j), 2), C)+term;
           
        end
        
    end
    
    
    if N == 2+i
        
        for j= 1:3
            
           
            term=tensorpower(paulimatrices(j),2)+term;
            
            
        end
        
    end
    
    if i>0
        
        term=kron(tensorpower(eye(2,2),i),term);
        
    end
    