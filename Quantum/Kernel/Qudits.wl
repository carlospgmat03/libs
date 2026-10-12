BeginPackage["Quantum`"]

(* {{{ Gate application *)
QuditApplyGate::usage = "QuditApplyGate[psi, U, targets, d] applies the d^k x d^k \
matrix U to the k qudits of local dimension d selected by the bitmask targets \
(qudit p = bit p, 0-indexed from the rightmost site). psi is a state vector. \
Qudit version of ApplyMultiQubitGate."

(* }}} Gate application *)
(* {{{ Reduced state *)
QuditReducedState::usage = "QuditReducedState[psi, leave, d] returns the reduced density \
matrix of the pure state psi on the qudits selected by the bitmask leave (qudit p = bit p, \
0-indexed from the rightmost site), tracing out the rest. Qudit version of \
PartialTrace[Proyector[psi], leave]."

(* }}} Reduced state *)
(* {{{ Tensor product *)
QuditTensorProduct::usage = "QuditTensorProduct[loc, psi1, psi2, d] returns the tensor \
product of the qudit state vectors psi1 and psi2 (local dimension d), with psi2 on the \
qudits selected by the bitmask loc and psi1 on the rest (qudit p = bit p, 0-indexed from \
the rightmost site). Qudit version of tensorProduct[loc, psi1, psi2]."

(* }}} Tensor product *)

Begin["`Private`"]

(* {{{ Gate application *)
QuditApplyGate::dimerr =
  "Dimensions of U (`1`) and psi (`2`) are inconsistent with targets `3` and d = `4`.";
(* {{{ *) quditOrder[mask_, n_] :=
  (* indice de arreglo j <-> posicion n - j; los sitios de mask van al frente *)
  Module[{idx = Sort[n - Flatten[Position[Reverse[IntegerDigits[mask, 2, n]], 1] - 1]]},
    Join[idx, Complement[Range[n], idx]]]
(* }}} *)
(* {{{ *) quditToFront[psi_, mask_, n_, d_] :=
  (* Devuelve la matriz d^k x d^(n-k) con los sitios de mask en las filas,
     y la permutacion de indices. *)
  Module[{k = DigitCount[mask, 2, 1], order = quditOrder[mask, n]},
    {ArrayReshape[Transpose[ArrayReshape[psi, ConstantArray[d, n]], InversePermutation[order]],
       {d^k, d^(n - k)}], order}]
(* }}} *)
(* {{{ *) QuditApplyGate[psi_?VectorQ, U_?MatrixQ, targets_Integer, d_Integer] :=
  Module[{n, mat, order},
    n = Log[d, Length[psi]];
    If[!IntegerQ[n] || Length[U] =!= d^DigitCount[targets, 2, 1] || targets >= 2^n,
      Message[QuditApplyGate::dimerr, Dimensions[U], Length[psi], targets, d];
      Return[$Failed]];
    {mat, order} = quditToFront[psi, targets, n, d];
    Flatten[Transpose[ArrayReshape[U . mat, ConstantArray[d, n]], order]]]
(* }}} *)
(* }}} Gate application *)
(* {{{ Reduced state *)
QuditReducedState::dimerr =
  "Length of psi (`1`) is inconsistent with mask `2` and d = `3`.";
(* {{{ *) QuditReducedState[psi_?VectorQ, leave_Integer, d_Integer] :=
  Module[{n, mat},
    n = Log[d, Length[psi]];
    If[!IntegerQ[n] || leave >= 2^n,
      Message[QuditReducedState::dimerr, Length[psi], leave, d]; Return[$Failed]];
    mat = First[quditToFront[psi, leave, n, d]];
    mat . ConjugateTranspose[mat]]
(* }}} *)
(* }}} Reduced state *)
(* {{{ Tensor product *)
QuditTensorProduct::dimerr =
  "Lengths `1`, `2` are inconsistent with loc `3` and d = `4`.";
(* {{{ *) QuditTensorProduct[loc_Integer, psi1_?VectorQ, psi2_?VectorQ, d_Integer] :=
  Module[{n1, n2, n},
    n1 = Log[d, Length[psi1]]; n2 = Log[d, Length[psi2]]; n = n1 + n2;
    If[!IntegerQ[n1] || !IntegerQ[n2] || DigitCount[loc, 2, 1] =!= n2 || loc >= 2^n,
      Message[QuditTensorProduct::dimerr, Length[psi1], Length[psi2], loc, d];
      Return[$Failed]];
    (* psi2 al frente, psi1 atras; luego se regresa a su lugar *)
    Flatten[Transpose[ArrayReshape[Flatten[KroneckerProduct[psi2, psi1]], ConstantArray[d, n]],
      quditOrder[loc, n]]]]
(* }}} *)
(* }}} Tensor product *)

End[]

EndPackage[]
