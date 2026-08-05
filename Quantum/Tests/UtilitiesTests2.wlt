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
