(* ::Package:: *)

(* {{{ *) BeginPackage["QuantumDavid`",{"Carlos`", "Quantum`"}]
Isingterm::usage = "Isingterm[i_,j_,N_]"
IsingChain::usage = "IsingChain[J_,N_]"
Hallvsall::usage = "Hallvsall[J_,N_]"
IsingChainInhom::usage = "IsingChainInhom[J_,Jinhom_,N_]"
sigma::usage = "sigma[i_,qubit_,N_] Operador de Pauli i aplicado al qubit con etiqueta qubit, para un sistema con un total de N qubits"
HK::usage = "HK[N_,bx_,bz_]"
matrixU::usage = "matrixU[bx_,qubits_,topology_]"
IPRSym::usage = "IPRSym[bx_,wk_,topology_]"
PRSym::usage = "PRSym[bx_,wk_,topology_]"
NumberToBinary::usage = "NumberToBinary[u_,bits_]"
ToBinary::usage = "ToBinary[state_]"
ToBase::usage = "ToBase[list_]"
K::usage = "K[qubits_]"
testsym::usage = "testsym[bx_,qubits_,steps_]"
Extractbyk::usage = "Extractbyk[k_,{values_,vecs_}]]"
IPRbyCohstateSymbetter::usage= "IPRbyCohstateSymbetter[\[Theta]_,\[Phi]_,list_,dim_]"
vecsk::usage= "vecsk[qubits_,k_]"
IPRbyCohstateSym::usage= "IPRbyCohstateSym[\[Theta]_, \[Phi]_, bx_, {values_, vecs_}, 
  topology_] "
ModifiedCoherentState::usage = "ModifiedCoherentState[\[Theta]_, \[Phi]_, qubits_]"
ModifiedCoherentState2::usage = "ModifiedCoherentState2[\[Theta]_, \[Phi]_, qubits_]"
IPRbyStatebetter::usage = "IPRbyStatebetter[stateinput_,list_,vecsk_]"
StateToDirac::usage = "VectorViewer[vec_] It shows the vector in Dirac notation in qubit representation."

Begin["Private`"] 

Isingterm[i_,j_,N_]:=Module[{list},
list=Table[If[k==i||k==j,PauliMatrix[3],IdentityMatrix[2]],{k,0,N-1}];
Apply[KroneckerProduct,list]
];

IsingChain[J_,N_]:=J*Sum[Isingterm[i,i+1,N],{i,0,N-2}]+J*Isingterm[N-1,0,N];

Hallvsall[J_,N_]:=Module[{i,j,H},
H=ConstantArray[0,{2^N,2^N}];
For[i=0,i<N,i++,
For[j=1+i,j<N,j++,
H=H+J*Isingterm[i,j,N];
]
];
H
];

IsingChainInhom[J_,Jinhom_,N_]:=J*Sum[Isingterm[i,i+1,N],{i,1,N-2}]+J*Isingterm[N-1,0,N]+Jinhom*Isingterm[0,1,N];

sigma[i_,qubit_,N_]:=Module[{list},
list=Table[If[k==qubit,PauliMatrix[i],IdentityMatrix[2]],{k,0,N-1}];
Apply[KroneckerProduct,list]
];

HK[N_,bx_,bz_]:=bx Sum[sigma[1,qubit,N],{qubit,0,N-1}]+bz Sum[sigma[3,qubit,N],{qubit,0,N-1}];

matrixU[bx_,qubits_,topology_]:=Module[{HKi,HI},
If[topology==4,HKi=HK[qubits,bx,1.4]+sigma[1,0,qubits]\[Delta]bx,HKi=HK[qubits,bx,1.4]];
Switch[topology,1,HI=IsingChain[1.0,qubits],2,HI=Hallvsall[1.0,qubits],3,HI=IsingChainInhom[1.0,1.0+\[Delta]J,qubits],4,HI=IsingChain[1.0,qubits]];
MatrixExp[-1.0*I HKi].MatrixExp[-1.0*I HI]
];

matrixU[bx_,qubits_,topology_,\[Delta]_]:=Module[{HKi,HI},
If[topology==4,HKi=HK[qubits,bx,1.4]+sigma[1,0,qubits]\[Delta],HKi=HK[qubits,bx,1.4]];
Switch[topology,1,HI=IsingChain[1.0,qubits],2,HI=Hallvsall[1.0,qubits],3,HI=IsingChainInhom[1.0,1.0+\[Delta],qubits],4,HI=IsingChain[1.0,qubits]];
MatrixExp[-1.0*I HKi].MatrixExp[-1.0*I HI]
];

IPRSym[bx_,wk_,topology_,\[Delta]_]:=Module[{U,list,qubits,U0},
qubits=Log[2,Length[Transpose[wk][[1]]]];
U=matrixU[bx,qubits,topology,\[Delta]];
U0=Dagger[wk].U.wk;
list=Orthogonalize[Eigenvectors[U0]];
1/Length[list]Total[Abs[list]^4,2]
];

PRSym[bx_,wk_,topology_,\[Delta]_]:=Module[{U,list,qubits,U0},
qubits=Log[2,Length[Transpose[wk][[1]]]];
U=matrixU[bx,qubits,topology,\[Delta]];
U0=Dagger[wk].U.wk;
list=Orthogonalize[Eigenvectors[U0]];
Total[Table[Total[Abs[list[[index]]]^4]^(-1),{index,1,Length[Transpose[wk]]}]]
];

NumberToBinary[u_,bits_]:=Module[{uu,out},uu=u;Reverse[Table[out=Mod[uu,2];uu=IntegerPart[uu/2];out,{bits}]]];

ToBinary[state_]:=NumberToBinary[Position[state,1][[1]][[1]]-1,Log[2,Length[state]]];

ToBase[list_]:=Module[{sum},
sum=1;
Table[If[list[[i]]==1,sum=sum+2^(Length[list]-i)],{i,Length[list]}];
SparseArray[sum->1,2^Length[list]]//Normal
];

K[qubits_]:=Module[{B},
B=ConstantArray[0,2^(2*qubits)];
Table[If[Mod[i,2^qubits+2]==0,B[[i+1]]=1,B[[i+1]]=0],{i,0,2^(2*qubits)/2-1}];
Table[If[Mod[i,2^qubits+2]==1,B[[i+2^(2*qubits)/2+1]]=1,B[[i+2^(2*qubits)/2+1]]=0],{i,0,2^(2*qubits)/2-1}];
Partition[B,2^qubits]
];

testsym[bx_,qubits_,steps_]:=Module[{A,U,sta,values,vecs,list},
U=matrixU[bx,qubits,1];
A=N[K[qubits]];
{values,vecs}=Eigensystem[A];
vecs=Orthogonalize[vecs];
list=Flatten[Table[sta=MatrixPower[U,steps].vecs[[i]];Chop[A.sta-values[[i]]sta],{i,1,2^qubits}]];
list=DeleteCases[list,0];
If[Length[list]==0,Print["No hay pex"],Print["Preocupate mi cabron"]]
];

Extractbyk[k_,{values_,vecs_}]:=Module[{pos,dim},
dim=Length[vecs[[1]]];
pos=Flatten[Position[Round[Chop[Log[values]*Log[2,dim]/(2*I Pi)]],k]];
Table[vecs[[i]],{i,pos}]
];

IPRbyCohstateSym[\[Theta]_, \[Phi]_, bx_, {values_, vecs_}, 
  topology_] := Module[{vecs0, w0, U0, U, sta, qubits, list, dim},
  qubits = Log[2, Length[vecs[[1]]]];
  U = matrixU[bx, qubits, topology];
  vecs0 = Extractbyk[0, {values, Orthogonalize[vecs]}];
  dim = vecs0 // Length;
  w0 = Transpose[vecs0];
  sta = Table[
    Chop[QuantumDotProduct[vecs0[[i]], 
      CoherentState[\[Theta], \[Phi], qubits]]], {i, dim}];
  U0 = Dagger[w0].U.w0;
  list = Orthogonalize[Eigenvectors[U0]];
  1/dim Total[
    Table[Abs[QuantumDotProduct[list[[i]], sta]]^4, {i, 1, dim}]]
  ];

IPRbyCohstateSymbetter[\[Theta]_,\[Phi]_,list_,vecsk_]:=Module[{state},
state=Table[Chop[QuantumDotProduct[vecsk[[i]],CoherentState[\[Theta],\[Phi],Log[2,Length[vecsk[[1]]]]]]],{i,Length[vecsk]}];
1/Length[vecsk] Total[Table[Abs[QuantumDotProduct[list[[i]],state]]^4,{i,1,Length[vecsk]}]]
];

IPRbyStatebetter[stateinput_,list_,vecsk_]:=Module[{state},
state=Table[Chop[QuantumDotProduct[vecsk[[i]],stateinput]],{i,Length[vecsk]}];
1/Length[vecsk] Total[Table[Abs[QuantumDotProduct[list[[i]],state]]^4,{i,1,Length[vecsk]}]]
];

vecsk[qubits_,k_]:=Module[{values,vecs},
{values, vecs} = Eigensystem[N[K[qubits]]];
Extractbyk[k, {values, Orthogonalize[vecs]}]
];

ModifiedCoherentState[\[Theta]_, \[Phi]_, qubits_] := 
 Flatten[KroneckerProduct[CoherentState[\[Theta], \[Phi], 3], 
   PauliMatrix[1].CoherentState[\[Theta], \[Phi], 1], 
   CoherentState[\[Theta], \[Phi], qubits - 4]], 1];

ModifiedCoherentState2[\[Theta]_, \[Phi]_, qubits_] := 
 Flatten[KroneckerProduct[CoherentState[\[Theta], \[Phi], 2], 
   KroneckerProduct[PauliMatrix[1].CoherentState[\[Theta], \[Phi], 1],
     PauliMatrix[1].CoherentState[\[Theta], \[Phi], 1]], 
   CoherentState[\[Theta], \[Phi], qubits - 4]], 1];

StateToDirac[vec_]:=Module[{vec2},
Quiet[
Chop[
Total[
DeleteCases[
Table[vec2=ConstantArray[0,Length[vec]];
vec2[[i]]=Round[Abs[Sign[vec[[i]]]]];
vec[[i]]*"|"<>ToString[TableForm[{ToBinary[vec2]}, TableSpacing->{1.2,1.2}],StandardForm]<>"\[RightAngleBracket]",{i,1,Length[vec]}],0]
]
]
]
];

End[] 
EndPackage[]
