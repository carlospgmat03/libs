BeginPackage["Quantum`"]

Purity::usage = "Purity[rho] returns the purity Tr[rho.rho] of a density matrix."

vonNeumannEntropy::usage = "vonNeumannEntropy[rho] returns the von Neumann entropy \
-Tr[rho Log2[rho]] of a density matrix. vonNeumannEntropy[eigenvalues] computes it \
from a list of eigenvalues."

Concurrence::usage = "Concurrence[rho] returns the concurrence of a 2-qubit density matrix. \
Concurrence[psi] for a pure 2-qubit state returns |<psi| sigma_y\[CircleTimes]sigma_y |psi*>|."

StateToBlochSphere::usage = "StateToBlochSphere[rho] returns the Bloch vector \
{Tr[rho.sigma_x], Tr[rho.sigma_y], Tr[rho.sigma_z]} for a single-qubit density matrix."

StateToBlochBasisRepresentation::usage = "StateToBlochBasisRepresentation[rho] returns the \
4x4 matrix T_{ij} = Tr[rho . sigma_i \[CircleTimes] sigma_j] for a 2-qubit density matrix \
(i,j = 0..3)."

QuantumMutualInformationReducedEntropies::usage = "QuantumMutualInformationReducedEntropies[rho] \
returns S(rho_A) + S(rho_B) - S(rho) for a 2-qubit density matrix."

QuantumMutualInformationMeasurements::usage = "QuantumMutualInformationMeasurements[rho] \
returns the measurement-based quantum mutual information (optimized over projective measurements \
on subsystem B) for a 2-qubit density matrix."

Discord::usage = "Discord[rho] returns the quantum discord of a 2-qubit density matrix, \
defined as QuantumMutualInformationReducedEntropies - QuantumMutualInformationMeasurements."

Begin["`Private`"]

(* {{{ Purity and entropy *)

(* {{{ *) Purity[rho_] := Tr[rho . rho]
(* }}} *)

(* {{{ *) vonNeumannEntropy[r_?MatrixQ] := vonNeumannEntropy[Eigenvalues[r]]
vonNeumannEntropy[evals_?VectorQ] := -Total[If[# == 0, 0, # Log[2, #]] & /@ Re[evals]]
(* }}} *)

(* }}} *)

(* {{{ Entanglement *)

sigmaysigmay = {{0,0,0,-1},{0,0,1,0},{0,1,0,0},{-1,0,0,0}};

(* {{{ *) Concurrence[rho_?MatrixQ] :=
  Module[{lambda},
    lambda = Sqrt[Abs[Eigenvalues[rho . sigmaysigmay . Conjugate[rho] . sigmaysigmay]]];
    Max[{2 Max[lambda] - Total[lambda], 0}]]

Concurrence[psi_?VectorQ] :=
  Abs[QuantumDotProduct[sigmaysigmay . Conjugate[psi], psi]]
(* }}} *)

(* }}} *)

(* {{{ Bloch sphere representations *)

(* {{{ *) StateToBlochSphere[rho_?MatrixQ] :=
  Chop[Tr /@ ({Pauli[1], Pauli[2], Pauli[3]} . rho)]
(* }}} *)

(* {{{ *) StateToBlochBasisRepresentation[rho_?MatrixQ] :=
  Table[Tr[rho . Pauli[{i, j}]], {i, 0, 3}, {j, 0, 3}]
(* }}} *)

(* Private helper for Discord/QuantumMutualInformationMeasurements *)
(* {{{ *) blochNormalFormVectors[rho_?MatrixQ] :=
  Module[{R, T, oa, cd, ob, A},
    R = Chop[StateToBlochBasisRepresentation[rho]];
    T = R[[2;;, 2;;]];
    {oa, cd, ob} = SingularValueDecomposition[T];
    A = Chop[DirectSum[{{1}}, Transpose[oa/Det[oa]]] . R . DirectSum[{{1}}, ob/Det[ob]]];
    {A[[2;;, 1]], A[[1, 2;;]], Diagonal[A][[2;;]]}]
(* }}} *)

(* }}} *)

(* {{{ Quantum mutual information and discord *)

(* {{{ *) QuantumMutualInformationReducedEntropies[rho_?MatrixQ] :=
  Module[{ra = PartialTrace[rho, 1], rb = PartialTrace[rho, 2]},
    vonNeumannEntropy[ra] + vonNeumannEntropy[rb] - vonNeumannEntropy[rho]]
(* }}} *)

(* {{{ *) QuantumMutualInformationMeasurements[rho_?MatrixQ] :=
  Module[{th, phi},
    Maximize[qmiMeasure[rho, {th, phi}], {th, phi}][[1]]]

qmiMeasure[rho_?MatrixQ, {th_, phi_}] :=
  vonNeumannEntropy[PartialTrace[rho, 2]] - conditionalEntropy[rho, {th, phi}]

conditionalEntropy[rho_, {th_, phi_}] :=
  Module[{x, a, b, c, mp, mm, ps, lambdas},
    x = {2 Cos[th] Sin[th] Cos[phi], 2 Cos[th] Sin[th] Sin[phi], 2 Cos[th]^2 - 1};
    {a, b, c} = blochNormalFormVectors[rho];
    mp = a + c x;  mm = a - c x;
    ps = {(1 + b . x)/2, (1 - b . x)/2};
    lambdas = {{(1 + Norm[mp]/(1 + b.x))/2, (1 - Norm[mp]/(1 + b.x))/2},
               {(1 + Norm[mm]/(1 - b.x))/2, (1 - Norm[mm]/(1 - b.x))/2}};
    ps[[1]] vonNeumannEntropy[lambdas[[1]]] + ps[[2]] vonNeumannEntropy[lambdas[[2]]]]
(* }}} *)

(* {{{ *) Discord[rho_?MatrixQ] :=
  QuantumMutualInformationReducedEntropies[rho] - QuantumMutualInformationMeasurements[rho]
(* }}} *)

(* }}} *)

End[]

EndPackage[]
