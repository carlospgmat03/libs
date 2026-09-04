(* Tests for new Utilities.wl additions *)

(* Dagger *)
VerificationTest[
  Dagger[{{1, I}, {0, 1}}],
  {{1, 0}, {-I, 1}},
  TestID -> "Dagger-complex-matrix"]

VerificationTest[
  Dagger[IdentityMatrix[3]],
  IdentityMatrix[3],
  TestID -> "Dagger-identity-is-self-adjoint"]

VerificationTest[
  Dagger[Dagger[{{1, 2+I}, {3, 4-I}}]],
  {{1, 2+I}, {3, 4-I}},
  TestID -> "Dagger-involution"]

(* Commutator *)
(* [sigma_x, sigma_y] = 2I sigma_z *)
VerificationTest[
  Normal[Commutator[Pauli[1], Pauli[2]]],
  Normal[2 I Pauli[3]],
  TestID -> "Commutator-PauliX-PauliY-gives-2i-PauliZ"]

VerificationTest[
  Commutator[IdentityMatrix[2], {{1,2},{3,4}}],
  {{0,0},{0,0}},
  TestID -> "Commutator-identity-commutes-with-anything"]

(* DirectSum *)
VerificationTest[
  DirectSum[{{1, 0}, {0, -1}}, {{2}}],
  {{1, 0, 0}, {0, -1, 0}, {0, 0, 2}},
  TestID -> "DirectSum-2x2-and-1x1"]

VerificationTest[
  DirectSum[{{1}}, {{1}}],
  {{1, 0}, {0, 1}},
  TestID -> "DirectSum-two-1x1-gives-2x2-identity"]

VerificationTest[
  DirectSum[{{{1,0},{0,1}}, {{2,0},{0,2}}, {{3}}}],
  {{1,0,0,0,0},{0,1,0,0,0},{0,0,2,0,0},{0,0,0,2,0},{0,0,0,0,3}},
  TestID -> "DirectSum-list-form-three-blocks"]

(* xlog2x *)
VerificationTest[xlog2x[0], 0, TestID -> "xlog2x-zero-gives-zero"]
VerificationTest[xlog2x[1], 0, TestID -> "xlog2x-one-gives-zero"]
VerificationTest[xlog2x[1/2], -(1/2), TestID -> "xlog2x-half"]

(* ValidDensityMatrix *)
VerificationTest[
  ValidDensityMatrix[{{1/2, 0}, {0, 1/2}}],
  True,
  TestID -> "ValidDensityMatrix-maximally-mixed-qubit-is-valid"]

VerificationTest[
  ValidDensityMatrix[{{1, 0}, {0, 0}}],
  True,
  TestID -> "ValidDensityMatrix-pure-state-is-valid"]

VerificationTest[
  ValidDensityMatrix[{{2, 0}, {0, 0}}],
  False,
  TestID -> "ValidDensityMatrix-trace-2-is-invalid"]

VerificationTest[
  ValidDensityMatrix[{{3/4, 1/2}, {1/2, 1/4}}],
  False,
  TestID -> "ValidDensityMatrix-negative-eigenvalue-is-invalid"]

(* GetDigitQuantumPosition *)
(* 5 = 101 in binary; bit 0 (LSB) = 1, bit 1 = 0, bit 2 = 1 *)
VerificationTest[GetDigitQuantumPosition[5, 0], 1, TestID -> "GetDigitQuantumPosition-5-bit0"]
VerificationTest[GetDigitQuantumPosition[5, 1], 0, TestID -> "GetDigitQuantumPosition-5-bit1"]
VerificationTest[GetDigitQuantumPosition[5, 2], 1, TestID -> "GetDigitQuantumPosition-5-bit2"]
VerificationTest[GetDigitQuantumPosition[5, 3], 0, TestID -> "GetDigitQuantumPosition-5-beyond"]

(* GellMann: SU(n) has n^2-1 generators *)
VerificationTest[
  Length[GellMann[2]], 3,
  TestID -> "GellMann-SU2-has-3-generators"]

VerificationTest[
  Length[GellMann[3]], 8,
  TestID -> "GellMann-SU3-has-8-generators"]

VerificationTest[
  Dimensions[GellMann[3][[1]]], {3, 3},
  TestID -> "GellMann-SU3-generators-are-3x3"]

(* GellMann generators are Hermitian *)
VerificationTest[
  AllTrue[GellMann[3], Norm[Normal[#] - Normal[ConjugateTranspose[#]], "Frobenius"] < 10^-10 &],
  True,
  TestID -> "GellMann-SU3-generators-are-hermitian"]

(* GellMann generators are traceless *)
VerificationTest[
  AllTrue[GellMann[3], Abs[Tr[#]] < 10^-10 &],
  True,
  TestID -> "GellMann-SU3-generators-are-traceless"]

(* PermutationMatrices *)
VerificationTest[
  Length[PermutationMatrices[3]], 6,
  TestID -> "PermutationMatrices-3-elements-gives-6-matrices"]

VerificationTest[
  Dimensions[PermutationMatrices[3][[1]]], {3, 3},
  TestID -> "PermutationMatrices-each-matrix-is-3x3"]

(* ---------------- GellMann ---------------- *)

VerificationTest[
  Normal /@ GellMann[2],
  PauliMatrix /@ {1, 2, 3},
  TestID -> "GellMann-n2-is-Pauli-matrices"
]

VerificationTest[
  Length[GellMann[3]],
  8,
  TestID -> "GellMann-n3-has-n2-minus-1-matrices"
]

(* Standard generalized Gell-Mann normalization: Tr[T_a.T_b] == 2 delta_ab *)
VerificationTest[
  Chop[Table[Tr[GellMann[3][[a]] . GellMann[3][[b]]], {a, 8}, {b, 8}]],
  2 IdentityMatrix[8],
  TestID -> "GellMann-n3-trace-orthonormality"
]

VerificationTest[
  Chop[Tr /@ GellMann[3]],
  ConstantArray[0, 8],
  TestID -> "GellMann-n3-traceless"
]

VerificationTest[
  Chop[Total[Flatten[Abs[Dagger[#] - #]]] & /@ GellMann[3]],
  ConstantArray[0, 8],
  TestID -> "GellMann-n3-hermitian"
]

VerificationTest[
  GellMann[0],
  $Failed,
  {GellMann::posint},
  TestID -> "GellMann-zero-raises-message"
]

VerificationTest[
  GellMann[-3],
  $Failed,
  {GellMann::posint},
  TestID -> "GellMann-negative-raises-message"
]

(* ApplyPermutation: 0-indexed push convention, consistent with ApplySwap *)

(* identity permutation leaves vector unchanged *)
VerificationTest[
  ApplyPermutation[{1, 2, 3, 4}, {0, 1}],
  {1, 2, 3, 4},
  TestID -> "ApplyPermutation-identity-vector"]

(* perm={1,0}: swap qubits 0 and 1; |01> -> |10> *)
VerificationTest[
  ApplyPermutation[{0, 1, 0, 0}, {1, 0}],
  {0, 0, 1, 0},
  TestID -> "ApplyPermutation-swap-2qubit-vector"]

(* consistency: ApplyPermutation swap matches ApplySwap *)
VerificationTest[
  ApplyPermutation[{0, 1, 0, 0}, {1, 0}],
  ApplySwap[{0, 1, 0, 0}, 0, 1],
  TestID -> "ApplyPermutation-matches-ApplySwap-0-1"]

(* perm={2,1,0}: qubit 0->2, 1->1, 2->0; |001> -> |100> (same as ApplySwap[..,0,2]) *)
VerificationTest[
  ApplyPermutation[{0, 1, 0, 0, 0, 0, 0, 0}, {2, 1, 0}],
  {0, 0, 0, 0, 1, 0, 0, 0},
  TestID -> "ApplyPermutation-reverse-3qubit-|001>->|100>"]

VerificationTest[
  ApplyPermutation[{0, 1, 0, 0, 0, 0, 0, 0}, {2, 1, 0}],
  ApplySwap[{0, 1, 0, 0, 0, 0, 0, 0}, 0, 2],
  TestID -> "ApplyPermutation-matches-ApplySwap-0-2"]

(* cyclic perm={1,2,0}: qubit 0->1, 1->2, 2->0; |010> -> |100> *)
VerificationTest[
  ApplyPermutation[{0, 0, 1, 0, 0, 0, 0, 0}, {1, 2, 0}],
  {0, 0, 0, 0, 1, 0, 0, 0},
  TestID -> "ApplyPermutation-cyclic-3qubit-|010>->|100>"]

(* inverse cyclic perm={2,0,1}: qubit 0->2, 1->0, 2->1; |010> -> |001> *)
VerificationTest[
  ApplyPermutation[{0, 0, 1, 0, 0, 0, 0, 0}, {2, 0, 1}],
  {0, 1, 0, 0, 0, 0, 0, 0},
  TestID -> "ApplyPermutation-inv-cyclic-3qubit-|010>->|001>"]

(* density matrix: |01><01| -> |10><10| under swap *)
VerificationTest[
  ApplyPermutation[{{0,0,0,0},{0,1,0,0},{0,0,0,0},{0,0,0,0}}, {1, 0}],
  {{0,0,0,0},{0,0,0,0},{0,0,1,0},{0,0,0,0}},
  TestID -> "ApplyPermutation-density-matrix-swap"]

(* density matrix result matches ApplySwap result *)
VerificationTest[
  ApplyPermutation[{{0,0,0,0},{0,1,0,0},{0,0,0,0},{0,0,0,0}}, {1, 0}],
  ApplySwap[{{0,0,0,0},{0,1,0,0},{0,0,0,0},{0,0,0,0}}, 0, 1],
  TestID -> "ApplyPermutation-density-matrix-matches-ApplySwap"]

(* PartialTrace invariance: a unitary acting only on the traced-out qubits cannot
   change the reduced state of the kept qubits, since
   Tr_B[(1_A x U_B) rho (1_A x U_B^dag)] = Tr_B[rho] by cyclicity on the B factor.
   All cases use 4 qubits (dim 16): `mask` selects the kept qubits, the gate targets
   are the complementary bits. Seeds are fixed so a failure is reproducible. *)

(* pure state, contiguous split: keep qubits {0,1} (mask=3), U on qubits {2,3} (=12) *)
VerificationTest[
  Module[{psi, u},
    SeedRandom[20260904];
    psi = RandomState[16]; u = CUEMember[4];
    Norm[PartialTrace[ApplyMultiQubitGate[psi, u, 12], 3]
         - PartialTrace[psi, 3], "Frobenius"] < 10^-10],
  True,
  TestID -> "PartialTrace-invariant-under-unitary-on-traced-qubits-vector"]

(* pure state, interleaved split: keep qubits {0,2} (mask=5), U on qubits {1,3} (=10).
   This exercises the MergeTwoIntegers/ExtractDigits bit interleaving. *)
VerificationTest[
  Module[{psi, u},
    SeedRandom[20260905];
    psi = RandomState[16]; u = CUEMember[4];
    Norm[PartialTrace[ApplyMultiQubitGate[psi, u, 10], 5]
         - PartialTrace[psi, 5], "Frobenius"] < 10^-10],
  True,
  TestID -> "PartialTrace-invariant-under-unitary-on-traced-qubits-vector-interleaved"]

(* mixed state: exercises the MatrixQ branch of PartialTrace, which is separate code *)
VerificationTest[
  Module[{rho, u},
    SeedRandom[20260906];
    rho = RandomDensityMatrix[16]; u = CUEMember[4];
    Norm[PartialTrace[ApplyMultiQubitGate[rho, u, 10], 5]
         - PartialTrace[rho, 5], "Frobenius"] < 10^-10],
  True,
  TestID -> "PartialTrace-invariant-under-unitary-on-traced-qubits-matrix"]

(* asymmetric split: keep qubit {0} (mask=1), U on qubits {1,2,3} (=14), an 8x8 unitary *)
VerificationTest[
  Module[{psi, u},
    SeedRandom[20260907];
    psi = RandomState[16]; u = CUEMember[8];
    Norm[PartialTrace[ApplyMultiQubitGate[psi, u, 14], 1]
         - PartialTrace[psi, 1], "Frobenius"] < 10^-10],
  True,
  TestID -> "PartialTrace-invariant-under-unitary-on-traced-qubits-1-vs-3"]

(* Negative control: a unitary on the KEPT qubits does change the reduced state.
   Without this, the tests above would still pass if PartialTrace ignored `mask`
   or ApplyMultiQubitGate were a no-op. *)
VerificationTest[
  Module[{psi, u},
    SeedRandom[20260908];
    psi = RandomState[16]; u = CUEMember[4];
    Norm[PartialTrace[ApplyMultiQubitGate[psi, u, 5], 5]
         - PartialTrace[psi, 5], "Frobenius"] > 10^-6],
  True,
  TestID -> "PartialTrace-changes-under-unitary-on-kept-qubits"]
