BeginPackage["Quantum`"]

ApplyOperator::usage = "ApplyOperator[E, rho] returns E . rho . Dagger[E]."

ApplyChannel::usage = "ApplyChannel[{E1, E2, ...}, rho] applies the quantum channel with \
Kraus operators {E_k} to rho: Sum_k E_k . rho . Dagger[E_k]."

ApplyKrausOperators::usage = "ApplyKrausOperators[ops, rho] is an alias for ApplyChannel."

KrausOperatorsForSuperUnitaryEvolution::usage = "KrausOperatorsForSuperUnitaryEvolution[psienv, U] \
returns the Kraus operators for the channel induced by unitary U on a system-environment composite \
when the environment starts in state psienv (Niel-Chuang p. 360)."

Reshuffle::usage = "Reshuffle[Phi] returns the reshuffled matrix (realignment) of Phi. \
Useful for the Jamiolkowski isomorphism."

JamiolkowskiStateToOperatorChoi::usage = "JamiolkowskiStateToOperatorChoi[rho] converts a \
Choi state to the corresponding Choi operator (multiplies by Sqrt[dim])."

JamiolkowskiOperatorChoiToState::usage = "JamiolkowskiOperatorChoiToState[O] converts a Choi \
operator to the Choi state (divides by Sqrt[dim])."

ToSuperOperatorSpaceComputationalBasis::usage = "ToSuperOperatorSpaceComputationalBasis[rho] \
flattens a density matrix into a vector in superoperator (Liouville) space."

FromSuperOperatorSpaceComputationalBasis::usage = "FromSuperOperatorSpaceComputationalBasis[v] \
reshapes a superoperator-space vector back to a matrix."

ToSuperOperatorSpacePauliBasis::usage = "ToSuperOperatorSpacePauliBasis[rho] expands rho in \
the normalized Pauli basis for an n-qubit system."

FromSuperOperatorSpacePauliBasis::usage = "FromSuperOperatorSpacePauliBasis[v] reconstructs a \
density matrix from its Pauli-basis expansion vector."

TransformationMatrixPauliBasisToComputationalBasis::usage = "TransformationMatrixPauliBasisToComputationalBasis[] \
returns the 4x4 change-of-basis matrix from the Pauli basis to the computational (Choi) basis \
for a single qubit."

AveragePurityChannelPauliBasis::usage = "AveragePurityChannelPauliBasis[Lambda] returns the \
average output purity (Haar-averaged over pure input states) for a channel represented in the \
Pauli basis by the matrix Lambda."

RandomTracePreservingMapChoiBasis::usage = "RandomTracePreservingMapChoiBasis[] returns a \
random trace-preserving map in the Choi basis for a single-qubit channel."

MakeQuantumChannel::usage = "MakeQuantumChannel[gate, steps, env] constructs the Pauli-basis \
superoperator for the channel obtained by evolving gate[state] for steps steps, tracing over env. \
Returns a list of 4x4 channel matrices (one per step)."

Begin["`Private`"]

(* {{{ Channel application *)

(* {{{ *) ApplyOperator[e_, rho_] := e . rho . ConjugateTranspose[e]
(* }}} *)

(* {{{ *) ApplyChannel[ops_List, rho_] := Total[ApplyOperator[#, rho] & /@ ops]
(* }}} *)

(* {{{ *) ApplyKrausOperators[ops_List, rho_] := ApplyChannel[ops, rho]
(* }}} *)

(* {{{ *) KrausOperatorsForSuperUnitaryEvolution[psienv_?VectorQ, u_?MatrixQ] :=
  Module[{nEnv = Length[psienv], nTotal = Length[u], nSub},
    nSub = nTotal / nEnv;
    Table[
      Table[
        Sum[u[[nSub k + i + 1, nSub m + j + 1]] psienv[[m + 1]], {m, 0, nEnv - 1}],
        {i, 0, nSub - 1}, {j, 0, nSub - 1}],
      {k, 0, nEnv - 1}]]
(* }}} *)

(* }}} *)

(* {{{ Jamiolkowski isomorphism *)

(* {{{ *) Reshuffle[phi_?MatrixQ] :=
  Module[{dim = Round[Sqrt[Length[phi]]]},
    Table[
      With[{mn = IntegerDigits[mn0, dim, 2], munu = IntegerDigits[munu0, dim, 2]},
        phi[[FromDigits[{mn[[1]], munu[[1]]}, dim] + 1,
             FromDigits[{mn[[2]], munu[[2]]}, dim] + 1]]],
      {mn0, 0, dim^2 - 1}, {munu0, 0, dim^2 - 1}]]
(* }}} *)

(* {{{ *) JamiolkowskiStateToOperatorChoi[rho_?MatrixQ] := Sqrt[Length[rho]] Reshuffle[rho]
(* }}} *)

(* {{{ *) JamiolkowskiOperatorChoiToState[op_?MatrixQ] := Reshuffle[op] / Sqrt[Length[op]]
(* }}} *)

(* }}} *)

(* {{{ Superoperator space *)

(* {{{ *) ToSuperOperatorSpaceComputationalBasis[rho_?MatrixQ] := Flatten[rho]
(* }}} *)

(* {{{ *) FromSuperOperatorSpaceComputationalBasis[v_?VectorQ] :=
  Partition[v, Round[Sqrt[Length[v]]]]
(* }}} *)

(* {{{ *) ToSuperOperatorSpacePauliBasis[rho_?MatrixQ] :=
  Module[{nq = Round[Log[2, Length[rho]]]},
    Chop[Table[
      Tr[rho . Pauli[IntegerDigits[i, 4, nq]] / 2.^(nq/2)],
      {i, 0, Length[rho]^2 - 1}]]]
(* }}} *)

(* {{{ *) FromSuperOperatorSpacePauliBasis[v_?VectorQ] :=
  Module[{nq = Round[Log[2, Length[v]] / 2]},
    Sum[v[[i+1]] Pauli[IntegerDigits[i, 4, nq]] / 2.^(nq/2), {i, 0, Length[v] - 1}]]
(* }}} *)

(* {{{ *) TransformationMatrixPauliBasisToComputationalBasis[] :=
  {{1/Sqrt[2], 0, 0, 1/Sqrt[2]},
   {0, 1/Sqrt[2], -I/Sqrt[2], 0},
   {0, 1/Sqrt[2], I/Sqrt[2], 0},
   {1/Sqrt[2], 0, 0, -1/Sqrt[2]}}
(* }}} *)

(* }}} *)

(* {{{ Random channels and averages *)

(* {{{ *) AveragePurityChannelPauliBasis[lambda_] :=
  Total[Power[lambda[[All, 1]], 2]] / 2 +
  Power[Norm[lambda[[2;;, 2;;]], "Frobenius"], 2] / 6
(* }}} *)

(* {{{ *) RandomTracePreservingMapChoiBasis[] :=
  Module[{psi},
    psi = Total[MapThread[tensorProduct,
      {TwoRandomOrhtogonalStates[8], TwoRandomOrhtogonalStates[2]}]] / Sqrt[2];
    Reshuffle[2 PartialTrace[Proyector[psi], 3]]]
(* }}} *)

(* {{{ *) MakeQuantumChannel[gate_, steps_Integer, env_] :=
  Module[{sigma, ch, Y, X, z, o},
    Y = EvolvGate[gate, steps, env, {-I, 1}];
    z = EvolvGate[gate, steps, env, {0, 1}];
    o = EvolvGate[gate, steps, env, {1, 0}];
    X = EvolvGate[gate, steps, env, {1, 1}];
    sigma[2] = Y - z - o;
    sigma[1] = X - z - o;
    sigma[0] = z + o;
    sigma[3] = o - z;
    Do[ch[i, j] = (1/2) Table[Tr[PauliMatrix[i] . sigma[j][[k]]], {k, steps}],
       {i, 0, 3}, {j, 0, 3}];
    Chop[Table[Table[ch[i,j][[k]], {i, 0, 3}, {j, 0, 3}], {k, steps}]]]
(* }}} *)

(* }}} *)

End[]

EndPackage[]
