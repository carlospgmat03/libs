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
(* }}} *)
Begin["`Private`"] 
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
(* }}} *)
End[] 
EndPackage[]

