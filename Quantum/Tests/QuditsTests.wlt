(* Tests for Qudits.wl *)

(* QuditApplyGate *)
(* d = 2: same as ApplyMultiQubitGate, for contiguous and non-contiguous masks *)
VerificationTest[
  Module[{psi = RandomState[2^4]},
    Table[
      With[{U = CUEMember[2^DigitCount[m, 2, 1]]},
        Max[Abs[QuditApplyGate[psi, U, m, 2] - ApplyMultiQubitGate[psi, U, m]]] < 10^-10],
      {m, {1, 8, 3, 5, 9, 6, 13, 15}}]],
  ConstantArray[True, 8],
  TestID -> "QuditApplyGate-d2-equals-ApplyMultiQubitGate"]

(* d = 3: preserves the norm (4 qudits, non-contiguous mask) *)
VerificationTest[
  Module[{psi = RandomState[3^4]},
    Abs[Norm[QuditApplyGate[psi, CUEMember[9], 9, 3]] - 1] < 10^-10],
  True,
  TestID -> "QuditApplyGate-d3-preserves-norm"]

(* d = 3: on the two leftmost of three qudits it equals (U ⊗ 1) . psi *)
VerificationTest[
  Module[{psi = RandomState[3^3], U = CUEMember[9]},
    Max[Abs[QuditApplyGate[psi, U, 6, 3] - KroneckerProduct[U, IdentityMatrix[3]] . psi]] < 10^-10],
  True,
  TestID -> "QuditApplyGate-d3-leftmost-equals-Kronecker"]

(* QuditReducedState *)
(* d = 2: same as PartialTrace[Proyector[psi], leave], for all masks of 4 qubits *)
VerificationTest[
  Module[{psi = RandomState[2^4]},
    Table[
      Max[Abs[QuditReducedState[psi, m, 2] - PartialTrace[Proyector[psi], m]]] < 10^-10,
      {m, 1, 15}]],
  ConstantArray[True, 15],
  TestID -> "QuditReducedState-d2-equals-PartialTrace"]

(* d = 3: trace 1 on a non-contiguous mask *)
VerificationTest[
  Abs[Tr[QuditReducedState[RandomState[3^4], 9, 3]] - 1] < 10^-10,
  True,
  TestID -> "QuditReducedState-d3-trace-one"]

(* d = 3: product state, the reduced state on the leftmost qudit is |a><a| *)
VerificationTest[
  Module[{a = RandomState[3], b = RandomState[9]},
    Max[Abs[QuditReducedState[Flatten[KroneckerProduct[a, b]], 4, 3] - Proyector[a]]] < 10^-10],
  True,
  TestID -> "QuditReducedState-d3-product-state"]

(* QuditTensorProduct *)
(* d = 2: same as tensorProduct[loc, psi1, psi2], for all 2-qubit locations in 4 qubits
   (includes 33-like split masks such as 9) *)
VerificationTest[
  Module[{psi1 = RandomState[4], psi2 = RandomState[4]},
    Table[
      Max[Abs[QuditTensorProduct[m, psi1, psi2, 2] - tensorProduct[m, psi1, psi2]]] < 10^-10,
      {m, {3, 5, 6, 9, 10, 12}}]],
  ConstantArray[True, 6],
  TestID -> "QuditTensorProduct-d2-equals-tensorProduct"]

(* d = 2: the qubit case used for psi0, loc = 33 on 6 qubits *)
VerificationTest[
  Module[{psi1 = RandomState[16], psi2 = RandomState[4]},
    Max[Abs[QuditTensorProduct[33, psi1, psi2, 2] - tensorProduct[33, psi1, psi2]]] < 10^-10],
  True,
  TestID -> "QuditTensorProduct-d2-loc33"]

(* d = 3: the reduced states on loc and on its complement are |psi2><psi2| and |psi1><psi1| *)
VerificationTest[
  Module[{psi1 = RandomState[9], psi2 = RandomState[9], psi},
    psi = QuditTensorProduct[9, psi1, psi2, 3];
    {Max[Abs[QuditReducedState[psi, 9, 3] - Proyector[psi2]]] < 10^-10,
     Max[Abs[QuditReducedState[psi, 6, 3] - Proyector[psi1]]] < 10^-10}],
  {True, True},
  TestID -> "QuditTensorProduct-d3-marginals"]
