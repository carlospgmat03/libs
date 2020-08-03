__precompile__() # Este comando es para que julia precompile el paquete

module quantum

using LinearAlgebra

export projector, basisstate, random_state, random_state_stat

#Generic Quantum Mechanics

"""
Projector from either one or two states
"""
function projector(state)
    return state*state'
end

"""
Projector from either one or two states
"""
function projector(state1,state2)
    return state1*state2'
end

"""
Basis state
"""
function basisstate(i,n)
    aux=zeros(ComplexF64,n)
    aux[i]=1
    return aux
end

"""
Creates a random state
"""
function random_state(dim::Int=2)
    v=randn(dim,1)+randn(dim,1)im
    v=v/norm(v)
    return vec(v)
end

"""
Creates a random state statistically normalized
"""
function random_state_stat(dim::Int=2)
    v=randn(dim,1)+randn(dim,1)im
    v=v/sqrt(2*dim)
    return vec(v)
end

"""
From digits like in mathematica
"""
function fromdigits(digits,base)
    sum([digits[k]*base^(k-1) for k=1:length(digits)])
end


"""
    apply_unitary!(psi, u, target_qubit)
This function applies efficiently an arbitrary unitary "u" to the target qubit, and modifies the state psi. 
"""
function apply_unitary!(psi, u, target_qubit)
    number_of_qubits = trailing_zeros(length(psi))
    end_outer_counter = 2^(number_of_qubits-target_qubit-1)-1
    for counter_left_bits = 0:end_outer_counter
        number_left=counter_left_bits<< (target_qubit+1)
        end_right_counter = number_left + (1<<target_qubit)-1
        for counter_right_bits = number_left:end_right_counter
            counter_right_bits_1 = counter_right_bits + (1<<target_qubit)
            psi[counter_right_bits+1], psi[counter_right_bits_1+1]=u*[psi[counter_right_bits+1], psi[counter_right_bits_1+1]]
        end
    end
end

#Apply permutations

#Swaps
function applyswap(state::Array{Complex{Float64},2},target1::Int64,target2::Int64)::Array{Complex{Float64},2}
    len=size(state)[1]
    aux=zeros(ComplexF64,len,len)
        for i in range(1,stop=len)
            for j in range(1,stop=len)
                digits1=digits(i-1,base=2,pad=Int(log2(len)))
                digits2=digits(j-1,base=2,pad=Int(log2(len)))
                digits1p=deepcopy(digits1)
                digits2p=deepcopy(digits2)
                digits1p[target1]=digits1[target2]
                digits1p[target2]=digits1[target1]
                digits2p[target1]=digits2[target2]
                digits2p[target2]=digits2[target1]
                aux[i,j]=state[fromdigits(digits1p,2)+1,fromdigits(digits2p,2)+1]
        end
    end
    aux
end

function applyswap(state::Array{Complex{Float64},1},target1::Int64,target2::Int64)::Array{Complex{Float64},2}
    len=length(state)
    aux=zeros(ComplexF64,len)
    for i in range(1,stop=len)
        digits1=digits(i-1,base=2,pad=Int(log2(len)))
        digits1p=deepcopy(digits1)
        digits1p[target1]=digits1[target2]
        digits1p[target2]=digits1[target1]
        aux[i]=state[fromdigits(digits1p,2)+1]
        end
    projector(aux)
end


function applyswappure(state::Array{Complex{Float64},1},target1::Int64,target2::Int64)::Array{Complex{Float64},1}
    len=length(state)
    aux=zeros(ComplexF64,len)
    for i in range(1,stop=len)
        digits1=digits(i-1,base=2,pad=Int(log2(len)))
        digits1p=deepcopy(digits1)
        digits1p[target1]=digits1[target2]
        digits1p[target2]=digits1[target1]
        aux[i]=state[fromdigits(digits1p,2)+1]
        end
    aux
end

#Everypermutation under construction

# Spin chains



"""
apply_ising!(psi, J, target_qubit_1, target_qubit_2)
"""
function apply_ising!(psi, J, target_qubit_1, target_qubit_2)
    expJ=exp(-im*J)
    expJc=conj(expJ)
    for i = 0: length(psi)-1
        if testbit(i,target_qubit_1) $ testbit(i,target_qubit_2)
            psi[i+1]*=expJc
        else
            psi[i+1]*=expJ
        end
    end
end 


"""
apply_kick!(psi, b, target_qubit)
"""
function apply_kick!(psi, b, target_qubit)
    phi=norm(b)
    if phi!=0
        b_normalized=b/phi
        sigma_n=sigmas[1]*b_normalized[1]+sigmas[2]*b_normalized[2]+sigmas[3]*b_normalized[3]
        u=eye(2)*cos(phi)-1.0im*sigma_n*sin(phi)
        apply_unitary!(psi, u, target_qubit)
    end
end 


end