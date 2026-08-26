BeginPackage["Quantum`"]

BitsOn::usage = "BitsOn[n] counts the number of 1-bits in the binary representation of n."

MergeTwoIntegers::usage = "MergeTwoIntegers[na, nb, ndigits] merges two integers based \
on the bit positions of ndigits: 1-bit positions receive the bits of na, 0-bit positions \
receive the bits of nb. Useful for indexing tensor products."

ExtractDigits::usage = "ExtractDigits[n, mask] splits the bits of n into two integers: \
bits at 0-positions of mask form the first output, bits at 1-positions form the second. \
Inverse of MergeTwoIntegers."

tensorProduct::usage = "tensorProduct[A, B] returns the Kronecker product of two matrices, \
or the flattened Kronecker product of two vectors. tensorProduct[loc, A, B] interleaves \
the subsystems according to the 1-bit positions of the integer loc."

tensorPower::usage = "tensorPower[A, n] returns the n-fold tensor power \
A \[CircleTimes] A \[CircleTimes] ... \[CircleTimes] A."

ApplyMultiQubitGate::usage = "ApplyMultiQubitGate[state, gate, targets] applies a \
multi-qubit gate to the qubits selected by the bitmask targets (qubit k = bit k, \
0-indexed from the LSB). state can be a vector or a density matrix."

ApplySwap::usage = "ApplySwap[state, j, k] swaps qubits j and k (0-indexed from the LSB) \
in a state vector or density matrix. ApplySwap[state, {j, k}] is a convenience form."

QuantumDotProduct::usage = "QuantumDotProduct[psi1, psi2] returns the inner product \
<psi1|psi2> = Conjugate[psi1] . psi2."

PartialTrace::usage = "PartialTrace[rho, mask] traces out the qubits whose bits are 0 in \
mask, keeping the qubits whose bits are 1. rho can be a density matrix or a pure state \
vector. Qubits are 0-indexed from the LSB."

PartialTranspose::usage = "PartialTranspose[A, mask] partially transposes the subsystem \
selected by the 1-bits of mask (0-indexed from the LSB)."

PartialTransposeFirst::usage = "PartialTransposeFirst[A] partially transposes the first \
subsystem of a bipartite system with equal subsystem dimensions."

PartialTransposeSecond::usage = "PartialTransposeSecond[A] partially transposes the second \
subsystem of a bipartite system with equal subsystem dimensions."

Dagger::usage = "Dagger[M] returns the conjugate transpose (Hermitian conjugate) of M."

Commutator::usage = "Commutator[A, B] returns the commutator A.B - B.A."

DirectSum::usage = "DirectSum[A1, A2] returns the block-diagonal direct sum A1 \[CirclePlus] A2. \
DirectSum[{A1, A2, ...}] folds the direct sum over a list."

xlog2x::usage = "xlog2x[x] returns x Log[2, x], with the convention xlog2x[0] = 0."

ValidDensityMatrix::usage = "ValidDensityMatrix[rho] returns True if rho is a valid density \
matrix (unit trace to 10^-13, non-negative eigenvalues after Chop)."

GetDigitQuantumPosition::usage = "GetDigitQuantumPosition[index, k] returns the k-th bit of \
index (0-indexed from the LSB). Equivalent to BitAnd[BitShiftRight[index, k], 1]."

GellMann::usage = "GellMann[n] gives the n^2-1 generalized Gell-Mann \
matrices for SU(n) (construction from \
https://blog.suchideas.com/2016/04/sun-gell-mann-matrices-in-mathematica/), \
normalized so that Tr[GellMann[n][[a]].GellMann[n][[b]]] == 2 \
KroneckerDelta[a,b]. For n=2 this gives the three Pauli matrices. The \
result does not include the identity; add \
Sqrt[2/n] IdentityMatrix[n] yourself if you need a complete Hermitian \
basis. Results are memoized per n."

BaseForHermitianMatrices::usage = "BaseForHermitianMatrices[j, n] returns the j-th element \
of an orthonormal basis for n x n Hermitian matrices: diagonal projectors for j <= n, then \
symmetric and antisymmetric off-diagonal pairs."

PermutationMatrices::usage = "PermutationMatrices[n] returns the list of all n! permutation \
matrices for permutations of n elements (uses the built-in PermutationMatrix)."

ApplyPermutation::usage = "ApplyPermutation[state, perm] applies a qubit permutation to state \
(vector or density matrix). perm is a 0-indexed list {p_0, ..., p_{nq-1}} using push convention: \
p_k is the output qubit that receives input qubit k (0-indexed from the LSB). \
ApplyPermutation[state, {1,0}] is equivalent to ApplySwap[state, 0, 1]."

Begin["`Private`"]

(* {{{ Bitwise utilities *)

(* {{{ *) BitsOn[n_Integer] := DigitCount[n, 2, 1]
(* }}} *)

(* {{{ *) MergeTwoIntegers[na_Integer, nb_Integer, ndigits_Integer] :=
  Module[{len, mask, result},
    len    = Max[DigitCount[ndigits, 2, 0], BitLength[nb]] + BitsOn[ndigits];
    mask   = PadRight[Reverse[IntegerDigits[ndigits, 2]], len];
    result = PadRight[{}, len];
    result[[Flatten[Position[mask, 1]]]] = Reverse[IntegerDigits[na, 2, BitsOn[ndigits]]];
    result[[Flatten[Position[mask, 0]]]] = Reverse[IntegerDigits[nb, 2, len - BitsOn[ndigits]]];
    FromDigits[Reverse[result], 2]]
(* }}} *)

(* {{{ *) ExtractDigits[n_Integer, mask_Integer] :=
  Module[{len, bits, mbits},
    len   = IntegerLength[n, 2];
    bits  = IntegerDigits[n, 2, len];
    mbits = IntegerDigits[mask, 2, len];
    {FromDigits[Pick[bits, mbits, 0], 2],
     FromDigits[Pick[bits, mbits, 1], 2]}]
(* }}} *)

(* {{{ *) GetDigitQuantumPosition[index_, k_] :=
  Module[{digits = Reverse[IntegerDigits[index, 2]]},
    If[k >= Length[digits], 0, digits[[k + 1]]]]
(* }}} *)

(* }}} *)

(* {{{ Tensor products *)

(* {{{ *) tensorProduct[A_?MatrixQ, B_?MatrixQ] := KroneckerProduct[A, B]
(* }}} *)

(* {{{ *) tensorProduct[A_?VectorQ, B_?VectorQ] := Flatten[KroneckerProduct[A, B]]
(* }}} *)

(* {{{ *) tensorProduct[loc_Integer, A_?MatrixQ, B_?MatrixQ] :=
  Module[{len, rIdx, cIdx},
    len = Length[A] Length[B];
    Normal[SparseArray[Flatten[Table[
      rIdx = 1 + ExtractDigits[r - 1, loc];
      cIdx = 1 + ExtractDigits[c - 1, loc];
      {r, c} -> A[[rIdx[[1]], cIdx[[1]]]] B[[rIdx[[2]], cIdx[[2]]]],
      {r, len}, {c, len}]]]]]
(* }}} *)

(* {{{ *) tensorProduct[loc_Integer, A_?VectorQ, B_?VectorQ] :=
  Module[{len, idx},
    len = Length[A] Length[B];
    Normal[SparseArray[Table[
      idx = 1 + ExtractDigits[i - 1, loc];
      i -> A[[idx[[1]]]] B[[idx[[2]]]],
      {i, len}]]]]
(* }}} *)

(* {{{ *) tensorPower[A_, n_Integer /; n >= 1] := Nest[KroneckerProduct[A, #] &, A, n - 1]
(* }}} *)

(* }}} *)

(* {{{ Gate application *)

ApplyMultiQubitGate::dimerr =
  "Gate dimension `1` is inconsistent with targets `2`; need Length[gate] == 2^BitsOn[targets].";

(* {{{ *) ApplyMultiQubitGate[state_?VectorQ, gate_?MatrixQ, targets_Integer] :=
  Module[{dGate, dRest, pos, statenew},
    dGate = Length[gate];
    dRest = Length[state] / dGate;
    If[2^BitsOn[targets] =!= dGate,
      Message[ApplyMultiQubitGate::dimerr, dGate, targets]; Return[$Failed]];
    statenew = state;
    Do[
      pos = Table[MergeTwoIntegers[ti, u, targets], {ti, 0, dGate - 1}] + 1;
      statenew[[pos]] = gate . state[[pos]],
      {u, 0, dRest - 1}];
    Chop[statenew]]
(* }}} *)

(* {{{ *) ApplyMultiQubitGate[state_?MatrixQ, gate_?MatrixQ, targets_Integer] :=
  Module[{dGate, dRest, pos, pos2, gateH, statenew},
    dGate = Length[gate];
    dRest = Length[state] / dGate;
    If[2^BitsOn[targets] =!= dGate,
      Message[ApplyMultiQubitGate::dimerr, dGate, targets]; Return[$Failed]];
    gateH    = ConjugateTranspose[gate];
    statenew = state;
    Do[
      pos  = Table[MergeTwoIntegers[ti, u,  targets], {ti, 0, dGate - 1}] + 1;
      pos2 = Table[MergeTwoIntegers[ti, u2, targets], {ti, 0, dGate - 1}] + 1;
      statenew[[pos, pos2]] = gate . state[[pos, pos2]] . gateH,
      {u, 0, dRest - 1}, {u2, 0, dRest - 1}];
    Chop[statenew]]
(* }}} *)

(* {{{ *) ApplySwap[state_, j_Integer, k_Integer] /; j === k := state
ApplySwap[state_, j_Integer, k_Integer] :=
  ApplyMultiQubitGate[state, {{1,0,0,0},{0,0,1,0},{0,1,0,0},{0,0,0,1}}, 2^j + 2^k]
ApplySwap[state_, {j_Integer, k_Integer}] := ApplySwap[state, j, k]
(* }}} *)

(* {{{ *) ApplyPermutation[state_?VectorQ, perm_List] :=
  Module[{nq = Length[perm], n = Length[state], idx},
    idx = Table[
      With[{b = Reverse[IntegerDigits[j-1, 2, nq]]},
        FromDigits[Reverse[b[[perm + 1]]], 2] + 1],
      {j, n}];
    state[[idx]]]

ApplyPermutation[state_?MatrixQ, perm_List] :=
  Module[{nq = Length[perm], n = Length[state], idx},
    idx = Table[
      With[{b = Reverse[IntegerDigits[j-1, 2, nq]]},
        FromDigits[Reverse[b[[perm + 1]]], 2] + 1],
      {j, n}];
    state[[idx, idx]]]
(* }}} *)

(* }}} *)

(* {{{ Inner product *)

(* {{{ *) QuantumDotProduct[psi1_?VectorQ, psi2_?VectorQ] := Conjugate[psi1] . psi2
(* }}} *)

(* }}} *)

(* {{{ Partial trace and transpose *)

(* {{{ *) PartialTrace[rho_?MatrixQ, mask_Integer] :=
  Module[{dKeep, dTrace},
    dKeep  = 2^BitsOn[mask];
    dTrace = Length[rho] / dKeep;
    Table[
      Sum[rho[[MergeTwoIntegers[b1-1, a-1, mask]+1,
               MergeTwoIntegers[b2-1, a-1, mask]+1]], {a, dTrace}],
      {b1, dKeep}, {b2, dKeep}]]
(* }}} *)

(* {{{ *) PartialTrace[psi_?VectorQ, mask_Integer] :=
  Module[{dKeep, mat, extr},
    dKeep = 2^BitsOn[mask];
    mat = SparseArray[
      Table[extr = ExtractDigits[i-1, mask]; {extr[[2]]+1, extr[[1]]+1} -> psi[[i]],
        {i, Length[psi]}],
      {dKeep, Length[psi]/dKeep}];
    Normal[mat . ConjugateTranspose[mat]]]
(* }}} *)

(* {{{ *) PartialTranspose[A_?MatrixQ, mask_Integer] :=
  Table[
    A[[BitAnd[ij-1, BitNot[mask]] + BitAnd[kl-1, mask] + 1,
       BitAnd[kl-1, BitNot[mask]] + BitAnd[ij-1, mask] + 1]],
    {ij, Length[A]}, {kl, Length[A]}]
(* }}} *)

(* {{{ *) PartialTransposeFirst[A_?MatrixQ] :=
  Module[{n = Round[Sqrt[Length[A]]]},
    Table[
      With[{row = IntegerDigits[r, n, 2], col = IntegerDigits[c, n, 2]},
        A[[FromDigits[{col[[1]], row[[2]]}, n]+1,
           FromDigits[{row[[1]], col[[2]]}, n]+1]]],
      {r, 0, n^2-1}, {c, 0, n^2-1}]]
(* }}} *)

(* {{{ *) PartialTransposeSecond[A_?MatrixQ] :=
  Module[{n = Round[Sqrt[Length[A]]]},
    Table[
      With[{row = IntegerDigits[r, n, 2], col = IntegerDigits[c, n, 2]},
        A[[FromDigits[{row[[1]], col[[2]]}, n]+1,
           FromDigits[{col[[1]], row[[2]]}, n]+1]]],
      {r, 0, n^2-1}, {c, 0, n^2-1}]]
(* }}} *)

(* }}} *)

(* {{{ Matrix utilities *)

(* {{{ *) Dagger[m_] := ConjugateTranspose[m]
(* }}} *)

(* {{{ *) Unprotect[Commutator]
Commutator[a_?MatrixQ, b_?MatrixQ] := a . b - b . a
Protect[Commutator]
(* }}} *)

(* {{{ *) DirectSum[a_, b_] :=
  Module[{da = Dimensions[a], db = Dimensions[b]},
    ArrayFlatten[{{a, ConstantArray[0, {da[[1]], db[[2]]}]},
                  {ConstantArray[0, {db[[1]], da[[2]]}], b}}]]
DirectSum[list_List] := Fold[DirectSum, list[[1]], Rest[list]]
(* }}} *)

(* {{{ *) xlog2x[x_] := If[x == 0, 0, x Log[2, x]]
(* }}} *)

(* {{{ *) ValidDensityMatrix[rho_?MatrixQ] :=
  Abs[Tr[rho] - 1] < 10^-13 && AllTrue[Chop[Re[Eigenvalues[rho]]], # >= 0 &]
(* }}} *)

(* }}} *)

(* {{{ Basis *)

(* {{{ GellMann *) 
GellMann::posint = "Argument `1` must be a positive integer.";
GellMann[n_Integer?Positive] :=
 GellMann[n] = Flatten[Table[
    (* Symmetric case *)
    SparseArray[{{j, k} -> 1, {k, j} -> 1}, {n, n}]
    , {k, 2, n}, {j, 1, k - 1}], 1]~
   Join~Flatten[Table[
     (* Antisymmetric case *)
     SparseArray[{{j, k} -> -I, {k, j} -> +I}, {n, n}]
     , {k, 2, n}, {j, 1, k - 1}], 1]~
   Join~Table[
    (* Diagonal case *)
    Sqrt[2/l/(l + 1)] SparseArray[
      Table[{j, j} -> 1, {j, 1, l}]~Join~{{l + 1, l + 1} -> -l}, {n, n}]
    , {l, 1, n - 1}];
GellMann[n_] := (Message[GellMann::posint, n]; $Failed)

(* }}} *)

(* {{{ TriangularRoot and Triangular *)
triangular[n_] := n (n + 1) / 2
triangularRoot[n_] := (-1 + Sqrt[8 n + 1]) / 2
(* }}} *)

(* {{{ *) BaseForHermitianMatrices[j_Integer, n_Integer] :=
  If[j <= n,
    SparseArray[{j,j} -> 1, {n,n}],
    Module[{k = Ceiling[triangularRoot[j - n]],
            r = j - n - triangular[Ceiling[triangularRoot[j - n]] - 1]},
      SparseArray[{{n - k, n - k + r}, {n - k + r, n - k}} -> 1/Sqrt[2], {n,n}]]]
(* }}} *)

(* {{{ *) PermutationMatrices[n_Integer] := PermutationMatrix /@ Permutations[Range[n]]
(* }}} *)

(* }}} *)

End[]

EndPackage[]
