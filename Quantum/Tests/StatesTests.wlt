(* Tests for States.wl *)

(* RandomGaussianComplex *)
VerificationTest[Head[RandomGaussianComplex[]], Complex,
  TestID -> "RandomGaussianComplex-returns-complex"]

VerificationTest[
  SeedRandom[1]; Chop[RandomGaussianComplex[0, 0]],
  0,
  TestID -> "RandomGaussianComplex-zero-sigma-returns-mean"]

(* RandomState *)
VerificationTest[
  Length[RandomState[4]], 4,
  TestID -> "RandomState-correct-dimension"]

VerificationTest[
  Abs[Norm[RandomState[8]] - 1] < 10^-10,
  True,
  TestID -> "RandomState-is-normalized"]

(* TwoRandomOrhtogonalStates *)
VerificationTest[
  Length[TwoRandomOrhtogonalStates[4]], 2,
  TestID -> "TwoRandomOrhtogonalStates-returns-two-states"]

VerificationTest[
  Module[{pair}, pair = TwoRandomOrhtogonalStates[6];
    Abs[Conjugate[pair[[1]]] . pair[[2]]] < 10^-10],
  True,
  TestID -> "TwoRandomOrhtogonalStates-are-orthogonal"]

VerificationTest[
  Module[{pair}, pair = TwoRandomOrhtogonalStates[6];
    Abs[Norm[pair[[1]]] - 1] < 10^-10 && Abs[Norm[pair[[2]]] - 1] < 10^-10],
  True,
  TestID -> "TwoRandomOrhtogonalStates-are-normalized"]

(* Proyector *)
VerificationTest[
  Proyector[{1, 0}],
  {{1, 0}, {0, 0}},
  TestID -> "Proyector-|0><0|"]

VerificationTest[
  Proyector[{0, 1}],
  {{0, 0}, {0, 1}},
  TestID -> "Proyector-|1><1|"]

VerificationTest[
  Proyector[{1, 0}, {0, 1}],
  {{0, 1}, {0, 0}},
  TestID -> "Proyector-|0><1|"]

(* GUEMember *)
VerificationTest[
  Dimensions[GUEMember[4]], {4, 4},
  TestID -> "GUEMember-correct-dimensions"]

VerificationTest[
  Module[{m}, m = GUEMember[4]; Norm[m - ConjugateTranspose[m], "Frobenius"] < 10^-10],
  True,
  TestID -> "GUEMember-is-hermitian"]

(* GOEMember *)
VerificationTest[
  Dimensions[GOEMember[4]], {4, 4},
  TestID -> "GOEMember-correct-dimensions"]

VerificationTest[
  Module[{m}, m = GOEMember[4]; Norm[m - Transpose[m], "Frobenius"] < 10^-10],
  True,
  TestID -> "GOEMember-is-real-symmetric"]

(* CUEMember *)
VerificationTest[
  Dimensions[CUEMember[4]], {4, 4},
  TestID -> "CUEMember-correct-dimensions"]

VerificationTest[
  Module[{u}, u = CUEMember[4];
    Norm[u . ConjugateTranspose[u] - IdentityMatrix[4], "Frobenius"] < 10^-10],
  True,
  TestID -> "CUEMember-is-unitary"]

(* PUEMember *)
VerificationTest[
  Dimensions[PUEMember[4]], {4, 4},
  TestID -> "PUEMember-correct-dimensions"]

VerificationTest[
  Module[{m}, m = PUEMember[4]; Norm[m - ConjugateTranspose[m], "Frobenius"] < 10^-10],
  True,
  TestID -> "PUEMember-is-hermitian"]

VerificationTest[
  Module[{evs}, evs = Eigenvalues[PUEMember[4]];
    And @@ (Re[#] >= -0.5 - 10^-10 && Re[#] <= 0.5 + 10^-10 & /@ evs)],
  True,
  TestID -> "PUEMember-eigenvalues-in-range"]

(* RandomDensityMatrix *)
VerificationTest[
  Dimensions[RandomDensityMatrix[4]], {4, 4},
  TestID -> "RandomDensityMatrix-correct-dimensions"]

VerificationTest[
  Abs[Tr[RandomDensityMatrix[4]] - 1] < 10^-10,
  True,
  TestID -> "RandomDensityMatrix-unit-trace"]

VerificationTest[
  Module[{rho}, rho = RandomDensityMatrix[4];
    Norm[rho - ConjugateTranspose[rho], "Frobenius"] < 10^-10],
  True,
  TestID -> "RandomDensityMatrix-is-hermitian"]

VerificationTest[
  Min[Re[Eigenvalues[RandomDensityMatrix[4]]]] >= -10^-10,
  True,
  TestID -> "RandomDensityMatrix-is-positive-semidefinite"]

(* RandomMixedState *)
VerificationTest[
  Dimensions[RandomMixedState[4, 3]], {4, 4},
  TestID -> "RandomMixedState-correct-dimensions"]

VerificationTest[
  Abs[Tr[RandomMixedState[4, 3]] - 1] < 10^-10,
  True,
  TestID -> "RandomMixedState-unit-trace"]

(* RandomHermitianMatrix *)
VerificationTest[
  Dimensions[RandomHermitianMatrix[2]], {4, 4},
  TestID -> "RandomHermitianMatrix-2qubits-is-4x4"]

VerificationTest[
  Module[{m}, m = RandomHermitianMatrix[2];
    Norm[m - ConjugateTranspose[m], "Frobenius"] < 10^-10],
  True,
  TestID -> "RandomHermitianMatrix-is-hermitian"]

(* BasisState *)
VerificationTest[BasisState[0, 4], {1,0,0,0}, TestID -> "BasisState-0th-of-4"]
VerificationTest[BasisState[2, 4], {0,0,1,0}, TestID -> "BasisState-2nd-of-4"]

(* Bell *)
VerificationTest[Bell[], {1,0,0,1}/Sqrt[2], TestID -> "Bell[]-is-|00>+|11>-over-sqrt2"]
VerificationTest[Bell[2], {1,0,0,1}/Sqrt[2], TestID -> "Bell[2]-matches-Bell[]"]
VerificationTest[Length[Bell[3]], 9, TestID -> "Bell[3]-has-9-components"]
VerificationTest[Abs[Norm[Bell[3]] - 1] < 10^-10, True, TestID -> "Bell[3]-is-normalized"]

(* GHZ *)
VerificationTest[
  GHZ[2], {1,0,0,1}/Sqrt[2],
  TestID -> "GHZ[2]-is-|00>+|11>-over-sqrt2"]

VerificationTest[
  GHZ[3], {1,0,0,0,0,0,0,1}/Sqrt[2],
  TestID -> "GHZ[3]-has-endpoints-only"]

(* Wstate *)
VerificationTest[
  Wstate[2], {0, 1, 1, 0}/Sqrt[2],
  TestID -> "Wstate[2]-is-|01>+|10>-over-sqrt2"]

VerificationTest[
  Abs[Norm[Wstate[3]] - 1] < 10^-10, True,
  TestID -> "Wstate[3]-is-normalized"]

(* CoherentState *)
VerificationTest[
  CoherentState[0, 0, 2], {1., 0., 0., 0.},
  TestID -> "CoherentState-theta=0-is-|00>"]

VerificationTest[
  Length[CoherentState[Pi/3, Pi/4, 3]], 8,
  TestID -> "CoherentState-3qubits-has-8-components"]

(* BlochSphereToState *)
VerificationTest[
  BlochSphereToState[{0, 0}], {1, 0},
  TestID -> "BlochSphereToState-north-pole-is-|0>"]

VerificationTest[
  Chop[BlochSphereToState[{Pi, 0}]], {0, 1},
  TestID -> "BlochSphereToState-south-pole-is-|1>"]

(* Cartesian north pole: phi indeterminate but state is still {1,0} *)
VerificationTest[
  BlochSphereToState[{0, 0, 1}], {1, 0},
  TestID -> "BlochSphereToState-Cartesian-north-pole"]

(* Cartesian equator at phi=0: {1,0,0} -> theta=Pi/2, phi=0 -> {1,1}/Sqrt[2] *)
VerificationTest[
  Chop[BlochSphereToState[{1, 0, 0}] - {1, 1}/Sqrt[2]], {0, 0},
  TestID -> "BlochSphereToState-Cartesian-equator-phi=0"]
