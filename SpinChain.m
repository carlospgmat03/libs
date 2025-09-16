(* ::Package:: *)

(* {{{ *) BeginPackage["SpinChain`",{"Carlos`", "Quantum`"}]
(* {{{ Primitives *)
ApplyMagnetickKick::usage = 
  "ApplyMagnetickKick[state_, b_, target_] applies a magnetic kick with field strength b to the specified target qubit in the quantum state. \
ApplyMagnetickKick[state_, b_] applies the kick to all qubits in the state.";

ApplyIsing::usage = 
  "ApplyIsing[state_, J_, i_, j_] applies an Ising interaction of strength J between qubits i and j in the quantum state.";
(* }}} *)
(* {{{ Full topologies*)

ApplyIsingChain::usage = 
  "ApplyIsingChain[state_, J_] applies Ising interactions of strength J along a periodic chain topology to the quantum state.";
ApplyChain::usage = 
  "ApplyChain[state_, J_, b_] applies a magnetic kick with field strength b to all qubits and then de Ising interaction in a chain topology.";
ApplyInverseChain::usage=
  "ApplyInverseChain[state_, J_, b_] applies the inverse operator as ApplyChain[state_, J_, b_]. ";
ApplyIsingTorus::usage=
  "ApplyIsingTorus[state_, J_] applies Ising interactions of strength J along a periodic torus topology to the quantum state.";
ApplyChainStar::usage="ApplyChainStar[state_, Jenv_,Jint_, b_] Se hace la topologia de la estrella con el magnetic kick. There is a Jenv of the chain, that acts like and environment for a single qubit that is interacting maybe with all others via Jint and returns the new quantum state."
ApplyIsingStar::usage = 
  "ApplyIsingStar[state_?VectorQ, Jenv_, Jint_] applies Ising interactions of strength J in a star topology and returns the new quantum state.";
ApplyCommonEnvironment::usage = 
  "ApplyCommonEnvironment[state_, b_, Jenv_, Jcoupling_] applies a magnetic kick with field strength b to the quantum state, then applies Ising interactions with environment coupling Jenv and edge coupling Jcoupling, returning the new quantum state. Se refiere a la topologia (a) del PRL de n-body Bell en PRA.";
ApplyDephasingChain::usage="ApplyDephasingChain[psi0_, Delta_, Jenv_, benv_, Jinteraction_] Se hace la topologia de una cadena que solo hace dephasing"
ApplyIsingStarEnvironment::usage = 
  "ApplyIsingStarEnvironment[state_, Jenv_] applies Ising interactions of strength Jenv in the environment, which is the chain part of the topology  and returns the new quantum state.";
ApplyMagnetickKickStarEnvironment::usage="The kick part in all but the first qubit. Think of all but the first as an enviroment tipically."
ApplyIsingStarInteractionQubitEnvironment::usage="Se hace el kick pero solo la interaccion con el  environment"
ApplyInomogeneousChain::usage="ApplyInomogeneousChain[state_?VectorQ, J_, J10_, b_] Se hace la cadena inhomogenea donde J10 indica la interaccion ising entre el primero y segundo qubit"
ApplyIsingAllVsAll::usage = 
  "ApplyIsingAllVsAll[state_, J_] applies Ising interactions of strength J between all pairs of qubits and returns the new quantum state.";
(* }}} *)
(* {{{ Explicit Matrices *)
IsingMatrix::usage="Get the matrix for the Ising Interaction Sigma_i Sigma_j, or for sum Sigma_i Sigma_j. In the first case, call as IsingMatrix[{IsingPosition1_Integer, IsingPosition2_Integer}, Qubits_Integer], and in the second,  IsingMatrix[IsingPositions_List, Qubits_Integer]"
SpinChainIsingMatrix::usage="The Ising matrix that has to come in the IsingMatrix routine, for a perdic spin chain. Call as SpinChainIsingMatrix[Qubits_]  "
SpinGridIsingMatrix::usage="The Ising matrix that has to come in the IsingMatrix routine, for a toric spin grid. Call as SpinGridIsingMatrix[{nx_, ny_}]"
MatrixPauliMagneticField::usage="Matrix corresponding to the hamiltonian b.sum sigma_j"
HamiltonianMagenitcChain::usage="Matrix Corresponding to the continuous Periodic ising spin chain with magnetic field"
HamiltonianMagenitcGrid::usage="Matrix Corresponding to the continuous toric ising grid chain with magnetic field"
ApplyMagnetickKickInhom::usage="ApplyMagnetickKickInhom[state_, b_, binhom_], Se aplica el magnetick kick donde binhom es el campo del spin 0"
(* }}} *)
(* }}} *)
Begin["Private`"] 
(* {{{ Primitives*)
(* {{{ *) ApplyMagnetickKick[state_?VectorQ, b_, Target_] := 
 Module[{RotationMatrix, statenew, i, pos},
  RotationMatrix = MatrixExp[-I MatrixPauliMagneticField[b,1]];
  statenew = state;
  For[i = 0, i < Length[state]/2, i++,
   pos = {MergeTwoIntegers[0, i, Power[2, Target]], 
      MergeTwoIntegers[1, i, Power[2, Target]]} + 1;
   statenew[[pos]] = RotationMatrix.statenew[[pos]];];
  statenew]
(* }}} *)
(* {{{ *) ApplyMagnetickKick[state_?VectorQ, b_] := 
 Module[{Qubits, statenew, QubitToAdress},
  Qubits = Log[2, Length[state]];
  statenew = state;
  For[QubitToAdress = 0, QubitToAdress < Qubits, QubitToAdress++, 
   (*Print["en el loop QubitToAdress="<>ToString[QubitToAdress]];*)
   statenew = ApplyMagnetickKick[statenew, b, QubitToAdress]];
  statenew]
(* }}} *)
(* {{{ *) ApplyMagnetickKickStarEnvironment[state_?VectorQ, b_] := 
 Module[{Qubits, statenew, QubitToAdress},
  Qubits = Log[2, Length[state]];
  statenew = state;
  For[QubitToAdress = 1, QubitToAdress < Qubits, QubitToAdress++, 
   (*Print["en el loop QubitToAdress="<>ToString[QubitToAdress]];*)
   statenew = ApplyMagnetickKick[statenew, b, QubitToAdress]];
  statenew]
(* }}} *)
(* {{{ *) ApplyIsing[state_?VectorQ, J_, Position1_, Position2_] := 
 Module[{scalar, i, statenew},
  scalar = Exp[-I J];
  statenew = state;
  For[i = 0, i < Length[state], i++,
   If[BitGet[i, Position1] == BitGet[i, Position2], 
    statenew[[i + 1]] = scalar statenew[[i + 1]], 
    statenew[[i + 1]] = Conjugate[scalar] statenew[[i + 1]]]];
  statenew]
(* }}} *)
(* }}} *)
(* {{{ Full topologies *)
(* {{{ *) ApplyIsingChain[state_?VectorQ, J_] := Module[{Qubits, statenew, QubitToAdress,q},
  Qubits = Log[2, Length[state]];
  statenew=state;
  For[q=0, q<Qubits-1, q++, 
    statenew = ApplyIsing[statenew, J, q, q+1];
  ];
  statenew = ApplyIsing[statenew, J, 0 , Qubits-1];
  statenew]
(* }}} *)
(* {{{ *) ApplyInverseChain[state_?VectorQ, J_, b_] := Module[{statenew},
 statenew = state;
 statenew = ApplyMagnetickKick[statenew, -b];
 statenew = ApplyIsingChain[statenew, -J];
 statenew ]
(* }}} *)
(* {{{ *) ApplyChain[state_?VectorQ, J_, b_] := Module[{statenew},
 statenew = state;
 statenew = ApplyIsingChain[statenew, J];
 statenew = ApplyMagnetickKick[statenew, b];
 statenew ]
(* }}} *)
(* {{{ *) ApplyCommonEnvironment[state_?VectorQ, b_, Jenv_, Jcoupling_] := 
 Module[{statenew, J, qubits}, qubits = Log[2, Length[state]];
  J = Table[Jenv, {qubits}];
  J[[1]] = 0;
  J[[2]] = Jcoupling;
  J[[-1]] = Jcoupling;
  ApplyIsing[ApplyMagnetickKick[state, b], J]]
(* }}} *)
(* {{{ *) ApplyDephasingChain[psi0_, Delta_, Jenv_, benv_, Jinteraction_] := Module[{statenew},
  statenew = psi0;
 (*U2 interno del env, las ising y la interaccion con el medio*)
 statenew = ApplyIsingStar[statenew, Jenv, Jinteraction];
 (* En el qubit solo se hace en direccion z, en el resto de la cadena
    donde sea, para que sea integrable/caotica. 
 *)
 statenew = ApplyMagnetickKickStarEnvironment[statenew, benv];
 statenew = ApplyMagnetickKick[statenew, {0, 0, Delta/2}, 0];
 statenew]
(* }}} *)
(* {{{ *) ApplyChainStar[state_?VectorQ, Jenv_,Jint_, b_] := Module[{statenew},
 statenew = state;
 statenew = ApplyIsingStar[statenew, Jenv, Jint];
 statenew = ApplyMagnetickKick[statenew, b];
 statenew ]
(* }}} *)
(* {{{ *) ApplyIsingStar[state_?VectorQ, Jenv_, Jint_] := Module[{Qubits, statenew, QubitToAdress, q},
  Qubits = Log[2, Length[state]];
If[IntegerQ[Qubits]==False,Print["Error: The state does not correspond to a integer number of qubits"];Abort[]];
  statenew=state;
  statenew = ApplyIsingStarEnvironment[statenew, Jenv];
  statenew = ApplyIsingStarInteractionQubitEnvironment[statenew, Jint];
  statenew
  ];
(* }}} *)
(* {{{ *) ApplyIsingStarEnvironment[state_?VectorQ, Jenv_] := Module[{Qubits, statenew, QubitToAdress, q},
  Qubits = Log[2, Length[state]];
If[IntegerQ[Qubits]==False,Print["Error: The state does not correspond to a integer number of qubits"];Abort[]];
  statenew=state;
  For[q=1, q<Qubits-1, q++, statenew = ApplyIsing[statenew, Jenv, q , q+1]; ];
  statenew = ApplyIsing[statenew, Jenv, Qubits-1 , 1]; 
  statenew
  ];
(* }}} *)
(* {{{ *) ApplyIsingStarInteractionQubitEnvironment[state_?VectorQ, Jint_] := Module[{Qubits, statenew, QubitToAdress, q},
  Qubits = Log[2, Length[state]];
If[IntegerQ[Qubits]==False,Print["Error: The state does not correspond to a integer number of qubits"];Abort[]];
  statenew=state;
  For[q=1, q<Qubits, q++, statenew = ApplyIsing[statenew, Jint, 0 , q]; ];
  statenew
  ];
(* }}} *)
(* {{{ *) ApplyInomogeneousChain[state_?VectorQ, J_, J10_, b_] := Module[{Qubits, statenew, QubitToAdress, q},
  Qubits = Log[2, Length[state]];
If[IntegerQ[Qubits]==False,Print["Error: The state does not correspond to a integer number of qubits"];Abort[]];
  statenew=state;
statenew=ApplyIsing[statenew, J10, 0, 1];
  For[q=1, q<Qubits-1, q++, 
	statenew = ApplyIsing[statenew, J, q , q + 1]; 
	];
  statenew = ApplyIsing[statenew, J, 0 , Qubits-1];
	statenew = ApplyMagnetickKick[statenew, b];
  statenew
  ];
(* }}} *)
(* {{{ *) ApplyMagnetickKickInhom[state_, b_, binhom_] := Module[{finalstate, Qubits, i},
 Qubits = Log[2, Length[state]];
  finalstate = ApplyMagnetickKick[state, binhom, 0];
  For[i = 1, i < Qubits, i++, 
   finalstate = ApplyMagnetickKick[finalstate, b, i]];
  finalstate]
(* }}} *)
(* {{{ *)ApplyIsingAllVsAll[state_,J_]:=Module[{statenew, Qubits, i, j},
 Qubits = Log[2, Length[state]];
  statenew = state;
For[i=0,i<Qubits,i++,
For[j=1+i, j<Qubits,j++,
statenew=ApplyIsing[statenew,J,i,j];
]];
statenew]
(* }}} *)
(* }}} *)
(* {{{ Explicit Matrices *)
IsingMatrix[{IsingPosition1_Integer, IsingPosition2_Integer}, Qubits_Integer] := Pauli[Table[ If[l == Qubits - IsingPosition1 || l == Qubits - IsingPosition2, 3, 0], {l, Qubits}]];
IsingMatrix[IsingPositions_List, Qubits_Integer] := Module[{IsingPosition}, Sum[IsingMatrix[IsingPosition, Qubits], {IsingPosition, IsingPositions}]]
SpinChainIsingMatrix[Qubits_] := Table[{Mod[i, Qubits], Mod[i + 1, Qubits]}, {i, 0, Qubits - 1}]

SpinGridIsingMatrix[{nx_, ny_}] := Join[Table[{i, i + 1 + If[Mod[i, nx] == nx - 1, -nx, 0]}, {i, 0, nx ny - 1}], Table[{i, Mod[i + nx, nx ny]}, {i, 0, nx ny - 1}]]

MatrixPauliMagneticField[MagneticField_, Qubits_] := MagneticField.{SumSigmaX[Qubits], SumSigmaY[Qubits], SumSigmaZ[Qubits]}
HamiltonianMagenitcChain[MagneticField_, J_, Qubits_] := MatrixPauliMagneticField[MagneticField, Qubits] + J IsingMatrix[SpinChainIsingMatrix[Qubits], Qubits]
HamiltonianMagenitcGrid[MagneticField_, J_, {nx_, ny_}] := MatrixPauliMagneticField[MagneticField, nx ny] + J IsingMatrix[SpinGridIsingMatrix[{nx, ny}]]
(* }}} *)
End[] 
EndPackage[]
