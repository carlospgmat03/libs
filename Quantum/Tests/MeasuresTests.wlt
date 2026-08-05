(* Tests for Measures.wl *)

(* Purity *)
VerificationTest[
  Purity[{{1/2, 0}, {0, 1/2}}],
  1/2,
  TestID -> "Purity-maximally-mixed-qubit"]

VerificationTest[
  Purity[{{1, 0}, {0, 0}}],
  1,
  TestID -> "Purity-pure-state"]

VerificationTest[
  Abs[Purity[RandomDensityMatrix[4]] - 1] <= 1,
  True,
  TestID -> "Purity-density-matrix-in-valid-range"]

(* vonNeumannEntropy *)
VerificationTest[
  vonNeumannEntropy[{{1/2, 0}, {0, 1/2}}],
  1,
  TestID -> "vonNeumannEntropy-maximally-mixed-qubit-is-1"]

VerificationTest[
  vonNeumannEntropy[{{1, 0}, {0, 0}}],
  0,
  TestID -> "vonNeumannEntropy-pure-state-is-0"]

(* Maximally mixed 2-qubit: entropy = 2 *)
VerificationTest[
  vonNeumannEntropy[IdentityMatrix[4] / 4],
  2,
  TestID -> "vonNeumannEntropy-maximally-mixed-2qubit-is-2"]

(* vonNeumannEntropy from eigenvalue list *)
VerificationTest[
  vonNeumannEntropy[{1/2, 1/2}],
  1,
  TestID -> "vonNeumannEntropy-from-eigenvalues"]

VerificationTest[
  vonNeumannEntropy[{1, 0}],
  0,
  TestID -> "vonNeumannEntropy-pure-eigenvalues"]

(* Concurrence *)
(* Bell state has concurrence 1 *)
VerificationTest[
  Abs[Concurrence[Proyector[{1,0,0,1}/Sqrt[2]]] - 1] < 10^-10,
  True,
  TestID -> "Concurrence-Bell-state-is-1"]

(* Separable state |00> has concurrence 0 *)
VerificationTest[
  Abs[Concurrence[Proyector[{1,0,0,0}]]] < 10^-10,
  True,
  TestID -> "Concurrence-product-state-is-0"]

(* Pure state form *)
VerificationTest[
  Abs[Concurrence[{1,0,0,1}/Sqrt[2]] - 1] < 10^-10,
  True,
  TestID -> "Concurrence-vector-Bell-state-is-1"]

(* StateToBlochSphere *)
(* |0><0| -> {0, 0, 1} *)
VerificationTest[
  StateToBlochSphere[{{1, 0}, {0, 0}}],
  {0, 0, 1},
  TestID -> "StateToBlochSphere-|0><0|-is-north-pole"]

(* |1><1| -> {0, 0, -1} *)
VerificationTest[
  StateToBlochSphere[{{0, 0}, {0, 1}}],
  {0, 0, -1},
  TestID -> "StateToBlochSphere-|1><1|-is-south-pole"]

(* Maximally mixed: Bloch vector = {0,0,0} *)
VerificationTest[
  Chop[StateToBlochSphere[{{1/2, 0}, {0, 1/2}}]],
  {0, 0, 0},
  TestID -> "StateToBlochSphere-maximally-mixed-is-origin"]

(* StateToBlochBasisRepresentation *)
VerificationTest[
  Dimensions[StateToBlochBasisRepresentation[IdentityMatrix[4]/4]],
  {4, 4},
  TestID -> "StateToBlochBasisRepresentation-is-4x4"]

(* For maximally mixed 2-qubit, only T_{00} = 1, rest 0 *)
VerificationTest[
  Chop[StateToBlochBasisRepresentation[IdentityMatrix[4]/4]],
  {{1,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},
  TestID -> "StateToBlochBasisRepresentation-maximally-mixed"]


(* QuantumMutualInformationReducedEntropies *)
(* For product state |00>, mutual information = 0 *)
VerificationTest[
  Abs[QuantumMutualInformationReducedEntropies[Proyector[{1,0,0,0}]]] < 10^-10,
  True,
  TestID -> "QuantumMutualInformationReducedEntropies-product-state-is-0"]

(* For Bell state, mutual information = 2 (both reduced states maximally mixed) *)
VerificationTest[
  Abs[QuantumMutualInformationReducedEntropies[Proyector[{1,0,0,1}/Sqrt[2]]] - 2] < 10^-10,
  True,
  TestID -> "QuantumMutualInformationReducedEntropies-Bell-is-2"]

