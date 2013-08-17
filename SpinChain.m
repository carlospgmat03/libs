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
EvolvGate::usage="EvolvGate[Gate_, steps_, env_, state_]... Evoluciona cualquier estado aplicando un numero steps_ de veces la compuerta Gate de  la forma Gate[#, otherparameters_] donde debe ponerse # en el lugar donde Gate toma el estado"
MakeQuantumChannel::usage="MakeQuantumChannel[Gate_, steps_, env_] Donde Gate va de la forma Gate[#, otherparameters_]  donde debe ponerse # en el lugar donde Gate toma el estado"
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
(* {{{ *) ApplyChainStar[state_?VectorQ, Jenv_,Jint_, b_] := Module[{statenew},
 statenew = state;
 statenew = ApplyIsingStar[statenew, Jenv, Jint];
 statenew = ApplyMagnetickKick[statenew, b];
 statenew ]
(* }}} *)
(* }}} *)
(* {{{ *) ApplyIsingStar[state_?VectorQ, Jenv_, Jint_] := Module[{Qubits, statenew, QubitToAdress, q},
  Qubits = Log[2, Length[state]];
If[IntegerQ[Qubits]==False,Print["Error: The state does not correspond to a integer number of qubits"];Abort[]];
  statenew=state;
  For[q=1, q<Qubits, q++, 
  statenew = ApplyIsing[statenew, Jint, 0 , q];
];
(* }}} *)
(* }}} *)
(* {{{ *)
EvolvGate[Gate_, steps_, env_, state_]:=
 Module[{statefinal, list, gate, j},
  statefinal = tensorProduct[env,state];
  gate[statefinal_] := First[Gate & /@ {statefinal}];
  list = Table[statefinal = gate[statefinal];
    PartialTrace[statefinal, 1], {j, steps}
    ];
  list
  ]
(*{{{*)
(*}}}*)
(*}}}*)
MakeQuantumChannel[Gate_, steps_, env_] := 
 Module[{\[Sigma], channel, Y, X, cero, uno},
  Y = EvolvGate[Gate, steps, env, {-I, 1}];
  cero = EvolvGate[Gate, steps, env, {0, 1}];
  uno = EvolvGate[Gate, steps, env, {1, 0}];
  X = EvolvGate[Gate, steps, env, {1, 1}];
  \[Sigma][2] = Y - cero - uno;
  \[Sigma][1] = X - cero - uno;
  \[Sigma][0] = cero + uno;
  \[Sigma][3] = uno - cero;
  Table[channel[i, j] = 
    1/2 Table[Tr[PauliMatrix[i].\[Sigma][j][[k]]], {k, steps}], {i, 0,
     3}, {j, 0, 3}];
  Chop[Table[
    Table[channel[i, j][[k]], {i, 0, 3}, {j, 0, 3}], {k, steps}]]]
End[] 
EndPackage[]

