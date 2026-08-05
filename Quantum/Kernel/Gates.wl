BeginPackage["Quantum`"]

Pauli::usage = "Pauli[i] returns the i-th Pauli matrix (i = 0,1,2,3) as a sparse array. \
Pauli[{i1, i2, ...}] returns the tensor product Pauli[i1] \[CircleTimes] Pauli[i2] \[CircleTimes] ..."

Paulib::usage = "Paulib[{b1, b2, b3}] returns b1 Pauli[1] + b2 Pauli[2] + b3 Pauli[3]."

Begin["`Private`"]

Pauli[i_Integer] := SparseArray[PauliMatrix[i]]
Pauli[indices_List] := KroneckerProduct @@ (Pauli /@ indices)

Paulib[b_List /; Length[b] == 3] := b . Table[Pauli[i], {i, 1, 3}]

End[]

EndPackage[]
