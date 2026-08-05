BeginPackage["Quantum`"]

Pauli::usage = "Pauli[i] returns the i-th Pauli matrix (i = 0,1,2,3) as a sparse array. \
Pauli[{i1, i2, ...}] returns the tensor product Pauli[i1] \[CircleTimes] Pauli[i2] \[CircleTimes] ..."

Paulib::usage = "Paulib[{b1, b2, b3}] returns b1 Pauli[1] + b2 Pauli[2] + b3 Pauli[3]."

Hadamard::usage = "Hadamard[] returns the 2x2 Hadamard gate {{1,1},{1,-1}}/Sqrt[2]."

LocalOneQubitOperator::usage = "LocalOneQubitOperator[n, k, a] places the 2x2 operator a \
on qubit k of an n-qubit system (1-indexed from the MSB: k=1 is the leftmost/most-significant \
qubit). Returns a 2^n x 2^n matrix."

ControlNotMatrix::usage = "ControlNotMatrix[control, target, n] returns the n-qubit CNOT \
matrix with qubits 0-indexed from the LSB."

SWAPMatrix::usage = "SWAPMatrix[n, {j, k}] returns the 2^n x 2^n matrix that swaps qubits j \
and k (1-indexed from the MSB). Returns the identity when j == k. DEPRECATED: prefer \
ApplySwap which uses 0-indexed from LSB."

SumSigmaX::usage = "SumSigmaX[n] returns the total Pauli-X operator Sum_{k=1}^{n} sigma_x^{(k)} \
for an n-qubit system. Result is cached."

SumSigmaY::usage = "SumSigmaY[n] returns the total Pauli-Y operator Sum_{k=1}^{n} sigma_y^{(k)} \
for an n-qubit system. Result is cached."

SumSigmaZ::usage = "SumSigmaZ[n] returns the total Pauli-Z operator Sum_{k=1}^{n} sigma_z^{(k)} \
for an n-qubit system. Result is cached."

ApplyGate::usage = "ApplyGate[U, state, targetQubit] applies the 2x2 unitary U to the qubit \
targetQubit (0-indexed from the LSB) of state. Legacy: prefer ApplyMultiQubitGate."

ApplyControlGate::usage = "ApplyControlGate[U, state, targetQubit, controlQubit] applies U to \
targetQubit conditioned on controlQubit being 1 (both 0-indexed from the LSB). Legacy."

EvolvGate::usage = "EvolvGate[gate, steps, env, state] evolves state by applying gate \
(a function of the combined env\[CircleTimes]state system) for steps iterations, returning \
a list of reduced density matrices (traced over env) at each step."

Begin["`Private`"]

(* {{{ Pauli matrices *)

(* {{{ *) Pauli[i_Integer] := SparseArray[PauliMatrix[i]]
Pauli[{i_}] := Pauli[i]
Pauli[indices_List] := KroneckerProduct @@ (Pauli /@ indices)
(* }}} *)

(* {{{ *) Paulib[b_List /; Length[b] == 3] := b . Table[Pauli[i], {i, 1, 3}]
(* }}} *)

(* }}} *)

(* {{{ Single-qubit gates *)

(* {{{ *) Hadamard[] := {{1,1},{1,-1}} / Sqrt[2]
(* }}} *)

(* }}} *)

(* {{{ Multi-qubit operators *)

(* {{{ *) LocalOneQubitOperator[n_Integer, k_Integer, a_] /; 1 <= k <= n && Dimensions[a] == {2,2} :=
  KroneckerProduct[
    IdentityMatrix[2^(k-1), SparseArray],
    a,
    IdentityMatrix[2^(n-k), SparseArray]]
(* }}} *)

(* {{{ *) ControlNotMatrix[control_Integer, target_Integer, nQubits_Integer] :=
  SparseArray[
    Table[
      {i, flipTarget[control, target, i-1, nQubits] + 1},
      {i, 2^nQubits}] -> 1]

flipTarget[control_Integer, target_Integer, num_Integer, n_Integer] :=
  Module[{bits = IntegerDigits[num, 2, n]},
    If[bits[[-control-1]] == 1,
      bits[[-target-1]] = Mod[bits[[-target-1]] + 1, 2]];
    FromDigits[bits, 2]]
(* }}} *)

(* {{{ *) SWAPMatrix[n_Integer, {j_Integer, k_Integer}] /; 1 <= j <= n && 1 <= k <= n && j =!= k :=
  Total[LocalOneQubitOperator[n, j, #] . LocalOneQubitOperator[n, k, #] & /@ Table[Pauli[i], {i, 0, 3}]] / 2
SWAPMatrix[n_Integer, {j_Integer, k_Integer}] /; 1 <= j <= n && 1 <= k <= n && j === k :=
  IdentityMatrix[2^n, SparseArray]
(* }}} *)

(* {{{ *) SumSigmaX[n_Integer] := SumSigmaX[n] =
  Table[If[DigitCount[BitXor[i-1, j-1], 2, 1] == 1, 1, 0], {i, 2^n}, {j, 2^n}]
(* }}} *)

(* {{{ *) SumSigmaY[n_Integer] := SumSigmaY[n] =
  Table[If[DigitCount[BitXor[i-1, j-1], 2, 1] == 1, If[i > j, I, -I], 0], {i, 2^n}, {j, 2^n}]
(* }}} *)

(* {{{ *) SumSigmaZ[n_Integer] := SumSigmaZ[n] =
  DiagonalMatrix[Table[n - 2 DigitCount[i, 2, 1], {i, 0, 2^n - 1}]]
(* }}} *)

(* }}} *)

(* {{{ Legacy gate application *)

(* {{{ *) ApplyGate[u_, state_, targetQubit_] :=
  Module[{out = state, dGate = Length[u], dRest = Length[state]/Length[u], ie},
    Do[
      ie = MergeTwoIntegers[#, i, 2^targetQubit] & /@ Range[0, dGate-1];
      out[[ie + 1]] = u . state[[ie + 1]],
      {i, 0, dRest - 1}];
    out]
(* }}} *)

(* {{{ *) ApplyControlGate[u_, state_, targetQubit_, controlQubit_] :=
  Module[{out = state, normCtrl, iwc, ie},
    normCtrl = If[controlQubit > targetQubit, controlQubit - 1, controlQubit];
    Do[
      iwc = MergeTwoIntegers[1, i, 2^normCtrl];
      ie = MergeTwoIntegers[#, iwc, 2^targetQubit] & /@ {0, 1};
      out[[ie + 1]] = u . state[[ie + 1]],
      {i, 0, Length[state]/4 - 1}];
    out]
(* }}} *)

(* {{{ *) EvolvGate[gate_, steps_Integer, env_, state_] :=
  Module[{s = tensorProduct[env, state]},
    Table[s = gate[s]; PartialTrace[s, 1], {steps}]]
(* }}} *)

(* }}} *)

End[]

EndPackage[]
