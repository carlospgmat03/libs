(* {{{ *) BeginPackage["SpinChain`",{"Carlos`", "Quantum`"}]
(* {{{ Primitives *)
ApplyMagnetickKick::usage = "ApplyMagnetickKick[state_, b_, Target_] or ApplyMagnetickKick[state_, b_]"
ApplyIsing::usage = "ApplyIsing[state_, J_, Position1_, Position2_]"
(* }}} *)
(* {{{ Full topologies*)
ApplyIsingTorus::usage="Se hace la topologia de un toro. Solo la parte de Ising"
ApplyIsingChain::usage="Se hace la topologia de una cadena. Solo la parte de Ising"
ApplyChain::usage="Se hace la topologia de una cadena."
ApplyInverseChain::usage="Se hace la topologia de una cadena pero hacia atras en el tiempo."
ApplyCommonEnvironment::usage="Se refiere a la topologia (a) del PRL de n-body Bell en PRA."
ApplyChainStar::usage="ApplyChainStar[state_, Jenv_,Jint_, b_] Se hace la topologia de la estrella con el magnetic kick"
ApplyIsingStar::usage="Se hace la topologia de una estrella, solo la parte de Ising"
(* }}} *)
(* {{{ Explicit Matrices *)
IsingMatrix::usage="Get the matrix for the Ising Interaction Sigma_i Sigma_j, or for sum Sigma_i Sigma_j. In the first case, call as IsingMatrix[{IsingPosition1_Integer, IsingPosition2_Integer}, Qubits_Integer], and in the second,  IsingMatrix[IsingPositions_List, Qubits_Integer]"
SpinChainIsingMatrix::usage="The Ising matrix that has to come in the IsingMatrix routine, for a perdic spin chain. Call as SpinChainIsingMatrix[Qubits_]  "
SpinGridIsingMatrix::usage="The Ising matrix that has to come in the IsingMatrix routine, for a toric spin grid. Call as SpinGridIsingMatrix[{nx_, ny_}]"
MatrixPauliMagneticField::usage="Matrix corresponding to the hamiltonian b.sum sigma_j"
HamiltonianMagenitcChain::usage="Matrix Corresponding to the continuous Periodic ising spin chain with magnetic field"
HamiltonianMagenitcGrid::usage="Matrix Corresponding to the continuous toric ising grid chain with magnetic field"
(* }}} *)
(* }}} *)
Begin["Private`"] 
(* {{{ Primitives*)
(* {{{ *) ApplyMagnetickKick[state_?VectorQ, b_, Target_] := 
 Module[{RotationMatrix, statenew, i, pos},
  RotationMatrix = MatrixExp[-I Paulib[b]];
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
(* {{{ *) ApplyIsingChain[state_?VectorQ, J_] := Module[{Qubits, statenew, QubitToAdress},
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
  For[q=1, q<Qubits, q++, 
  statenew = ApplyIsing[statenew, Jint, 0 , q];
]];
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

