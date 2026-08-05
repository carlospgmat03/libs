(* Tests for Channels.wl *)

(* ApplyOperator *)
(* X |0><0| X = |1><1| *)
VerificationTest[
  ApplyOperator[{{0,1},{1,0}}, {{1,0},{0,0}}],
  {{0,0},{0,1}},
  TestID -> "ApplyOperator-X-flips-|0><0|"]

(* Unitary: U rho U† preserves trace *)
VerificationTest[
  Abs[Tr[ApplyOperator[CUEMember[4], RandomDensityMatrix[4]]] - 1] < 10^-10,
  True,
  TestID -> "ApplyOperator-preserves-trace"]

(* ApplyChannel / ApplyKrausOperators *)
(* Depolarizing channel: identity Kraus operators preserve state *)
VerificationTest[
  ApplyChannel[{IdentityMatrix[2]}, {{1/2, 0}, {0, 1/2}}],
  {{1/2, 0}, {0, 1/2}},
  TestID -> "ApplyChannel-identity-Kraus"]

(* Completely depolarizing: 4 Kraus = Pauli/2 maps any state to I/2 *)
VerificationTest[
  Chop[ApplyChannel[Table[Pauli[i]/Sqrt[4], {i, 0, 3}], {{1,0},{0,0}}] - {{1/2,0},{0,1/2}}],
  {{0,0},{0,0}},
  TestID -> "ApplyChannel-completely-depolarizing"]

VerificationTest[
  ApplyKrausOperators[{IdentityMatrix[2]}, {{1,0},{0,0}}],
  ApplyChannel[{IdentityMatrix[2]}, {{1,0},{0,0}}],
  TestID -> "ApplyKrausOperators-same-as-ApplyChannel"]

(* Reshuffle *)
VerificationTest[
  Dimensions[Reshuffle[Proyector[{1,0,0,1}/Sqrt[2]]]],
  {4, 4},
  TestID -> "Reshuffle-output-is-4x4"]

(* Reshuffle is an involution: R(R(M)) = M *)
VerificationTest[
  Module[{rho = RandomDensityMatrix[4]},
    Chop[Reshuffle[Reshuffle[rho]] - rho]],
  Table[0, {4}, {4}],
  TestID -> "Reshuffle-is-involution"]

(* Jamiolkowski round-trip *)
VerificationTest[
  Module[{rho = RandomDensityMatrix[4]},
    Chop[JamiolkowskiOperatorChoiToState[JamiolkowskiStateToOperatorChoi[rho]] - rho]],
  Table[0, {4}, {4}],
  TestID -> "Jamiolkowski-round-trip"]

(* ToSuperOperatorSpaceComputationalBasis / From *)
VerificationTest[
  Module[{rho = RandomDensityMatrix[4]},
    Chop[FromSuperOperatorSpaceComputationalBasis[
           ToSuperOperatorSpaceComputationalBasis[rho]] - rho]],
  Table[0, {4}, {4}],
  TestID -> "SuperOperatorComputationalBasis-round-trip"]

(* ToSuperOperatorSpacePauliBasis: maximally mixed has only zero-th component nonzero *)
VerificationTest[
  Chop[Rest[ToSuperOperatorSpacePauliBasis[IdentityMatrix[2]/2]]],
  {0, 0, 0},
  TestID -> "ToSuperOperatorSpacePauliBasis-maximally-mixed-only-first-component"]

(* From/ToSuperOperatorSpacePauliBasis round-trip *)
VerificationTest[
  Module[{rho = RandomDensityMatrix[2]},
    Chop[Normal[FromSuperOperatorSpacePauliBasis[
           ToSuperOperatorSpacePauliBasis[rho]]] - rho]],
  {{0, 0}, {0, 0}},
  TestID -> "SuperOperatorPauliBasis-round-trip"]

(* TransformationMatrixPauliBasisToComputationalBasis *)
VerificationTest[
  Dimensions[TransformationMatrixPauliBasisToComputationalBasis[]],
  {4, 4},
  TestID -> "TransformationMatrix-is-4x4"]

(* The matrix should be unitary *)
VerificationTest[
  Module[{m = TransformationMatrixPauliBasisToComputationalBasis[]},
    Chop[m . ConjugateTranspose[m] - IdentityMatrix[4]]],
  Table[0, {4}, {4}],
  TestID -> "TransformationMatrix-is-unitary"]

(* KrausOperatorsForSuperUnitaryEvolution *)
(* Unitary evolution: single Kraus operator = U (projected on env state) *)
VerificationTest[
  Module[{u = CUEMember[4], psienv = {1, 0}, kraus},
    kraus = KrausOperatorsForSuperUnitaryEvolution[psienv, u];
    (* For pure env initial state, sum_k E_k†E_k = I *)
    Chop[Sum[ConjugateTranspose[kraus[[k]]] . kraus[[k]], {k, Length[kraus]}] -
         IdentityMatrix[2]]],
  {{0, 0}, {0, 0}},
  TestID -> "KrausOperatorsForSuperUnitaryEvolution-completeness"]
