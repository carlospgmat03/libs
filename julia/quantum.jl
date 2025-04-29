__precompile__() # Este comando es para que julia precompile el paquete


module quantum

using LinearAlgebra

export projector, basisstate, random_state, random_state_stat, base_state, fromdigits, apply_unitary!, applyswap, applyswappure, apply_ising!, apply_kick!, testbit
export sigma_x, sigma_y, sigma_z, sigmas, merge_two_integers, pauli, parity_operator, apply_multiqubit_gate, apply_multiqubit_gate!, state_to_dirac, partial_trace, base_2, original_integer, extract_digits

#Generic Quantum Mechanics

@doc "Projector from either one or two states"
function projector(state)
    return state*state'
end

@doc "Projector from either one or two states"
function projector(state1,state2)
    return state1*state2'
end

@doc "Old one still in use for basis states desde el 1"
function basisstate(i::Int,n::Int)
    aux=zeros(ComplexF64,n)
    aux[i]=1
    return aux
end

@doc "Construye los elementos de la base comenzando desde el 0"
function base_state(i::Int,dim::Int)
    psi=zeros(Complex{Float64},dim)
    psi[i+1]=1;
    return psi
end

"""
Creates a random state
"""
function random_state(dim::Int=2)::Vector{Complex{<:AbstractFloat}}
    v=randn(dim,1)+randn(dim,1)im
    v=v/norm(v)
    return vec(v)
end

@doc "Creates a random state statistically normalized"
function random_state_stat(dim::Int=2)
    v=randn(dim,1)+randn(dim,1)im
    v=v/sqrt(2*dim)
    return vec(v)
end

@doc "From digits like in mathematica"
function fromdigits(digits,base)
    sum([digits[k]*base^(k-1) for k=1:length(digits)])
end


@doc "apply_unitary!(psi, u, target_qubit) This function applies efficiently an arbitrary unitary 'u' to the target qubit, and modifies the state psi."
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
@doc "Old swap function, use apply_swap instead"
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

@doc "Old swap function, use apply_swap instead"
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

@doc "Old swap function, use apply_swap instead"
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

@doc "apply_swap(state::Vector{T}, target1::Int, target2::Int) This function applies a swap operation to a given state."
function apply_swap(state::Vector{T}, target1::Int, target2::Int) where T
    new_state = copy(state)
    particles = log2(length(state))|>Int
    for i in 0:length(state)-1
        digits1=base_2(i, pad=particles)
        digits1p=copy(digits1)
        digits1p[target1]=digits1[target2]
        digits1p[target2]=digits1[target1]
        new_state[i+1] = state[original_integer(digits1p)+1]
    end
    return new_state
end

@doc "apply_swap!(state::Vector{T}, target1::Int, target2::Int) This function applies a swap operation to a given state and modifies the state in place."
function apply_swap!(state::Vector{T}, target1::Int, target2::Int) where T
    particles = log2(length(state))|>Int
    for i in 0:length(state)-1
        digits1=base_2(i, pad=particles)
        digits1p=copy(digits1)
        digits1p[target1]=digits1[target2]
        digits1p[target2]=digits1[target1]
        state[i+1] = state[original_integer(digits1p)+1]
    end
end

@doc "apply_swap(state::Matrix{T}, target1::Int, target2::Int) This function applies a swap operation to a given state."
function apply_swap(state::Matrix{T}, target1::Int, target2::Int) where T
    new_state = copy(state)  # Shallow copy is sufficient
    particles = log2(size(state, 1)) |> Int
    for (i, j) in Base.Iterators.product(0:size(state, 1)-1, 0:size(state, 2)-1)
        digits1 = base_2(i, pad=particles)
        digits2 = base_2(j, pad=particles)
        
        # Swap target1 and target2 in the binary representation
        digits1p = copy(digits1)
        digits2p = copy(digits2)
        digits1p[target1] = digits1[target2]
        digits1p[target2] = digits1[target1]
        digits2p[target1] = digits2[target2]
        digits2p[target2] = digits2[target1]
        
        # Update the new_state matrix
        new_state[i+1, j+1] = state[original_integer(digits1p)+1, original_integer(digits2p)+1]
    end
    return new_state
end

@doc "apply_swap!(state::Matrix{T}, target1::Int, target2::Int) This function applies a swap operation to a given state and modifies the state in place."
function apply_swap!(state::Matrix{T}, target1::Int, target2::Int) where T
    particles = log2(size(state, 1)) |> Int
    for (i, j) in Base.Iterators.product(0:size(state, 1)-1, 0:size(state, 2)-1)
        digits1 = base_2(i, pad=particles)
        digits2 = base_2(j, pad=particles)
        
        # Swap target1 and target2 in the binary representation
        digits1p = copy(digits1)
        digits2p = copy(digits2)
        digits1p[target1] = digits1[target2]
        digits1p[target2] = digits1[target1]
        digits2p[target1] = digits2[target2]
        digits2p[target2] = digits2[target1]
        
        # Update the new_state matrix
        state[i+1, j+1] = state[original_integer(digits1p)+1, original_integer(digits2p)+1]
    end
end

#Everypermutation under construction

# Spin chains

sigma_x=[0. 1.; 1. 0.]
sigma_y=[0. -im; im 0]
sigma_z=[1. 0.;0. -1.]
identity = [1. 0.; 0 1.]
sigmas = Dict(0 => identity, 1 => sigma_x, 2 => sigma_y, 3 => sigma_z)

@doc "function testbit(n, bit) This function test if a a given bit in position  'bit' of number 'n' is on."
function testbit(n, bit)
    ~(n&(1<<bit)==0)
end 

@doc "apply_ising!(psi, J, target_qubit_1, target_qubit_2)"
function apply_ising!(psi, J, target_qubit_1, target_qubit_2)
    expJ=exp(-im*J)
    expJc=conj(expJ)
    for i = 0: length(psi)-1
        if testbit(i,target_qubit_1) ⊻ testbit(i,target_qubit_2)
            psi[i+1]*=expJc
        else
            psi[i+1]*=expJ
        end
    end
end 



@doc "apply_kick!(psi, b, target_qubit)"
function apply_kick!(psi, b, target_qubit)
    phi=norm(b)
    if phi!=0
        b_normalized=b/phi
        sigma_n=sigmas[1]*b_normalized[1]+sigmas[2]*b_normalized[2]+sigmas[3]*b_normalized[3]
        u=[[1,0] [0,1]]*cos(phi)-1.0im*sigma_n*sin(phi)
        apply_unitary!(psi, u, target_qubit)
    end
end 


@doc "merge_two_integers, same functioning as usual"
function merge_two_integers(a::Int, b::Int, mask::Int)::Int
    result = 0
    bit_position = 0

    while mask != 0 || a != 0 || b != 0
        if mask & 1 == 1
            result |= (a & 1) << bit_position
            a >>= 1
        else
            result |= (b & 1) << bit_position
            b >>= 1
        end
        mask >>= 1
        bit_position += 1
    end

    return result
end

@doc "pauli(index, target, particles) quite self-explanatory, index goes form 0 to 3, where 0 is identity"
function pauli(index, target, particles)
    list = []
    for i in 0:particles-1
        if i == target
            push!(list, sigmas[index])
        else
            push!(list, sigmas[0])
        end
    end
    return kron(list...)
end

@doc "Parity operator"
function parity_operator(particles)
    list = [sigmas[1] for _ in 1:particles]
    return kron(list...)
end

function original_integer(list)
    return parse(Int, join(list); base=2)
end

function base_2(integer; pad= nothing)
    if pad == nothing
        return reverse(digits(integer, base = 2))
    else
        return reverse(digits(integer, base = 2, pad = pad))
    end
end

@doc "apply_multiqubit_gate(state::Vector{T}, gate, target) This function applies a multiqubit gate to a given state."
function apply_multiqubit_gate(state::Vector{T}, gate, target) where T
    new_state = copy(state)
    dim_target = size(gate)[1]
    dim_untouched = length(state)/dim_target|>Int
    for index_untouched in 0:dim_untouched-1
        pos = [merge_two_integers(index_target, index_untouched, target)+1 for index_target in 0:dim_target-1]
        new_state[pos] = gate*new_state[pos]
    end
    return new_state
end

@doc "apply_multiqubit_gate!(state::Vector{T}, gate, target) This function applies a multiqubit gate to a given state and modifies the state in place."
function apply_multiqubit_gate!(state::Vector{T}, gate, target) where T
    dim_target = size(gate)[1]
    dim_untouched = length(state)/dim_target|>Int
    for index_untouched in 0:dim_untouched-1
        pos = [merge_two_integers(index_target, index_untouched, target)+1 for index_target in 0:dim_target-1]
        state[pos] = gate*state[pos]
    end
end

@doc "apply_multiqubit_gate(state::Matrix{T}, gate, target) This function applies a multiqubit gate to a given state."
function apply_multiqubit_gate(state::Matrix{T}, gate, target) where T
    dim_target = size(gate)[1]
    dim_untouched = size(state)[1]/dim_target|>Int
    new_state = copy(state)
    for (index_untouched1,index_untouched2) in Base.product(0:dim_untouched-1, 0:dim_untouched-1)
        pos1 = [merge_two_integers(index_target, index_untouched1, target)+1 for index_target in 0:dim_target-1]
        pos2 = [merge_two_integers(index_target, index_untouched2, target)+1 for index_target in 0:dim_target-1]
        new_state[pos1, pos2] = gate*new_state[pos1, pos2]*gate'
    end
    return new_state
end

@doc "apply_multiqubit_gate!(state::Matrix{T}, gate, target) This function applies a multiqubit gate to a given state and modifies the state in place."
function apply_multiqubit_gate!(state::Matrix{T}, gate, target) where T
    dim_target = size(gate)[1]
    dim_untouched = size(state)[1]/dim_target|>Int
    for (index_untouched1,index_untouched2) in Base.product(0:dim_untouched-1, 0:dim_untouched-1)
        pos1 = [merge_two_integers(index_target, index_untouched1, target)+1 for index_target in 0:dim_target-1]
        pos2 = [merge_two_integers(index_target, index_untouched2, target)+1 for index_target in 0:dim_target-1]
        state[pos1, pos2] = gate*state[pos1, pos2]*gate'
    end
end

@doc "state_to_dirac(vector::Vector{<:Number}) Converts a vector representing a quantum state into Dirac notation."
function state_to_dirac(vector::Vector{<:Number})
    # Get the indices and values of non-zero elements
    indices = findall(x -> x != 0, vector)
    if isempty(indices)
        error("The vector does not represent a valid computational basis state.")
    end

    # Convert each index to binary representation
    num_qubits = Int(log2(length(vector)))
    binary_representations = [reverse(digits(index - 1, base=2, pad=num_qubits)) for index in indices]

    # Construct the ket notation with coefficients
    terms = [(string(vector[index]) * "|" * join(binary) * "⟩") for (index, binary) in zip(indices, binary_representations)]
    ket = join(terms, " + ")
    return ket
end

@doc "state_to_dirac(matrix::Matrix{<:Number}) Converts a matrix representing a quantum state into Dirac notation."
function state_to_dirac(matrix::Matrix{<:Number})
    # Get the size of the matrix
    dim = size(matrix, 1)
    if dim != size(matrix, 2)
        error("The matrix must be square.")
    end

    # Determine the number of qubits
    num_qubits = Int(log2(dim))
    if 2^num_qubits != dim
        error("The matrix size must be a power of 2.")
    end

    # Construct the ket-bra representation
    terms = []
    for i in 1:dim
        for j in 1:dim
            if matrix[i, j] != 0
                ket = "|" * join(reverse(digits(i - 1, base=2, pad=num_qubits))) * "⟩"
                bra = "⟨" * join(reverse(digits(j - 1, base=2, pad=num_qubits))) * "|"
                push!(terms, string(matrix[i, j], ket, bra))
            end
        end
    end

    # Join terms with " + " for the final representation
    ketbra = join(terms, " + ")
    return ketbra
end

@doc "partial_trace(state::Matrix{T}, target) This function calculates the partial trace of a quantum state over a specified target qubit."
function partial_trace(state::Matrix{T}, target) where T
    # Get the size of the state
    dim_total = size(state, 1)

    dim_target = 2^sum(base_2(target))

    # Calculate the dimensions of the reduced state
    dim_trace = dim_total/dim_target|>Int
    reduced_state = zeros(T, dim_target, dim_target)

    # Perform the partial trace over the target qubit
    for (i, j) in Base.product(0:dim_target-1, 0:dim_target-1)
        for k in 0:dim_trace-1
            # Calculate the indices for the original state
            index1 = merge_two_integers(i, k, target) + 1
            index2 = merge_two_integers(j, k, target) + 1
            reduced_state[i+1, j+1] += state[index1, index2]
        end
    end
    return reduced_state
end

@doc "partial_trace(state::Vector{T}, target) This function calculates the partial trace of a quantum state over a specified target qubit."
function partial_trace(state::Vector{T}, target) where T
    # Get the size of the state
    dim_total = length(state)

    # Calculate the dimensions of the reduced state
    dim_target = 2^sum(base_2(target))
    dim_trace = dim_total/dim_target|>Int
    reduced_state = zeros(T, dim_target, dim_target)

    # Perform the partial trace over the target qubit
    for (i, j) in Base.product(0:dim_target-1, 0:dim_target-1)
        for k in 0:dim_trace-1
            # Calculate the indices for the original state
            index1 = merge_two_integers(i, k, target) + 1
            index2 = merge_two_integers(j, k, target) + 1
            reduced_state[i+1, j+1] += state[index1] * conj(state[index2])
        end
    end
    return reduced_state
end

@doc "extract_digits(in::Int, target::Int) This function extracts the digits of a number in a given base. Same as implemented by Carlos in Mathematica"
function extract_digits(in::Int, target::Int)
    inlist = base_2(in)
    targetlist = base_2(target, pad=length(inlist))
    println("inlist: ", inlist)
    println("targetlist: ", targetlist)
    indices0 = findall(x -> x == 0, targetlist)
    indices1 = findall(x -> x == 1, targetlist)
    return (original_integer(inlist[indices0]),quantum.original_integer(inlist[indices1]))
end


end