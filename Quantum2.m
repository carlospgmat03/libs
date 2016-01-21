(* ::Package:: *)

(* {{{ *) BeginPackage["Quantum2`",{"Carlos`", "Quantum`"}]
Isingterm::usage = "Isingterm[i_,j_,N_]"
IsingChain::usage = "IsingChain[J_,N_]"
Hallvsall::usage = "Hallvsall[J_,N_]"
IsingChainInhom::usage = "IsingChainInhom[J_,Jinhom_,N_]"
sigma::usage = "sigma[i_,qubit_,N_] Operador de Pauli i aplicado al qubit con etiqueta qubit, para un sistema con un total de N qubits"
HK::usage = "HK[N_,bx_,bz_]"
matrixU::usage = "matrixU[bx_,bz_,qubits_,topology_]"
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
CharlieMeasure::usage = "CharlieMeasure[list_] or CharlieMeasure[list_]"
CharlieMeasureAve::usage = "CharlieMeasureAve[list_]"
CharlieMeasureForShowThings::usage = "CharlieMeasureForShowThings[list_]"
StairCase::usage = "StairCase[x_,eigen_]"
NNS::usage = "NNS[eigen_]"
Unfold::usage = "Unfold[list_]"
PBrody::usage = "PBrody[s_,q_]"
LSI::usage = "LSI[unfoldedlist_,bins_] Level Statistic Indicator or simply Gamma parameter"
G::usage = "Amplitude Damping Quantum Channel G[t_,\[Lambda]_,\[Omega]0_,\[Gamma]_]"
H::usage = "Binary Shannon Entropy H[p_]"
QuantumCapacityDamping::usage = "Quantum Capacity of the quantum damping, It must be specified the parameters of G, QuantumCapacityDamping[t_]"
ClassicalCapacityDamping::usage = "EA Classical Capacity of the quantum damping, It must be specified the parameters of G, QuantumCapacityDamping[t_]"
StepDecomposition::usage = "StepDecomposition[list_,\[Epsilon]_,elemnts_]"
BasisElement::usage = "BasisElement[i_,j_] Dont worry about this, not yet"
BasisElementOneIndex::usage = "BasisElementOneIndex[i_] Dont worry about this, not yet"
DivisivilityKindOf::usage = "DivisivilityKindOf[\[Lambda]1_,\[Lambda]2_,\[Lambda]3_] The lambdas state for the singular values of the unital channel
up rotations), then this function gives 0 when the channel is not CPTP, 1 if the channel is CPTP, 2 if it is p-divisible, 3 if it is compatible
with CP-divisible dynamics and 4 if the channel can be written as exp(L) with L Lindblad."


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

(*matrixU[bx_,qubits_,topology_]:=Module[{HKi,HI},
If[topology==4,HKi=HK[qubits,bx,1.4]+sigma[1,0,qubits]\[Delta]bx,HKi=HK[qubits,bx,1.4]];
Switch[topology,1,HI=IsingChain[1.0,qubits],2,HI=Hallvsall[1.0,qubits],3,HI=IsingChainInhom[1.0,1.0+\[Delta]J,qubits],4,HI=IsingChain[1.0,qubits]];
MatrixExp[-1.0*I HKi].MatrixExp[-1.0*I HI]
];*)

matrixU[bx_,bz_,qubits_,topology_]:=Module[{HKi,HI},
If[topology==4,HKi=HK[qubits,bx,bz]+sigma[1,0,qubits]\[Delta]bx,HKi=HK[qubits,bx,bz]];
Switch[topology,1,HI=IsingChain[1.0,qubits],2,HI=Hallvsall[1.0,qubits],3,HI=IsingChainInhom[1.0,1.0+\[Delta]J,qubits],4,HI=IsingChain[1.0,qubits]];
MatrixExp[-1.0*I HKi].MatrixExp[-1.0*I HI]
];

(*matrixU[bx_,qubits_,topology_,\[Delta]_]:=Module[{HKi,HI},
If[topology==4,HKi=HK[qubits,bx,1.4]+sigma[1,0,qubits]\[Delta],HKi=HK[qubits,bx,1.4]];
Switch[topology,1,HI=IsingChain[1.0,qubits],2,HI=Hallvsall[1.0,qubits],3,HI=IsingChainInhom[1.0,1.0+\[Delta],qubits],4,HI=IsingChain[1.0,qubits]];
MatrixExp[-1.0*I HKi].MatrixExp[-1.0*I HI]
];*)

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

Extractbyk[k_,{values_,vecs_}]:=Module[{pos,dim,logdim},
logdim=Log[2,dim];
dim=Length[vecs[[1]]];
pos=Flatten[Position[Round[Chop[Arg[values]*logdim/(2 Pi)]],k]];
Table[#[[i]],{i,pos}]&[vecs]
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
Total[Table[Abs[QuantumDotProduct[list[[i]],state]]^4,{i,1,Length[vecsk]}]]
];

IPRbyStatebetter[stateinput_,list_,vecsk_]:=Module[{state},
state=Table[Chop[QuantumDotProduct[vecsk[[i]],stateinput]],{i,Length[vecsk]}];
Total[Table[Abs[QuantumDotProduct[list[[i]],state]]^4,{i,1,Length[vecsk]}]]
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

CharlieMeasure[list2_] := 
  Module[{l, listD, Criticallistmin, Criticallistmax, position, maxi, 
    len, list},
   list = list2/Max[list2];
   l = Length[list];
   listD = Table[list[[i + 1]] - list[[i]], {i, l - 1}];
   Criticallistmax = 
    Sort[DeleteCases[
      Table[If[(listD[[i]] > 0 && 
           listD[[i + 1]] <= 0) || (listD[[i]] >= 0 && 
           listD[[i + 1]] < 0), {i, list[[i + 1]]}, 0], {i, l - 2}], 
      0], #1[[2]] > #2[[2]] &];
len=Length[Criticallistmax];
If[len==0,0,Max[{0,Max[Table[Criticallistmax[[i]][[2]]-Min[Take[list,Criticallistmax[[i]][[1]]]],{i,1,len}]]}]*Max[list2]]
];

CharlieMeasureAve[list2_] := 
  Module[{l, listD, Criticallistmin, Criticallistmax, position, maxi, 
    len, list},
   list = list2/Max[list2];
   l = Length[list];
   listD = Table[list[[i + 1]] - list[[i]], {i, l - 1}];
   Criticallistmax = 
    Sort[DeleteCases[
      Table[If[(listD[[i]] > 0 && 
           listD[[i + 1]] <= 0) || (listD[[i]] >= 0 && 
           listD[[i + 1]] < 0), {i, list[[i + 1]]}, 0], {i, l - 2}], 
      0], #1[[2]] > #2[[2]] &];
len=Length[Criticallistmax];
If[len==0,0,Max[{0,Max[Table[Criticallistmax[[i]][[2]]-Mean[Take[list,Criticallistmax[[i]][[1]]]],{i,1,len}]]}]*Max[list2]]
]

CharlieMeasureForShowThings[list_,deep_]:=Module[{l,listD,Criticallistmin,Criticallistmax,position,maxi,len,tab,positionmax,miningraph,maxingraph,minlist,positionmin},
l=Length[list];
listD=Table[list[[i+1]]-list[[i]],{i,l-1}];
Criticallistmax=Sort[DeleteCases[Table[If[(listD[[i]]>0&&listD[[i+1]]<=0)||(listD[[i]]>=0&&listD[[i+1]]<0),{i,list[[i+1]]},0],{i,l-2}],0],#1[[2]]>#2[[2]]&];
Criticallistmin=Table[If[(listD[[i]]<=0&&listD[[i+1]]>0)||(listD[[i]]<0&&listD[[i+1]]>=0),list[[i+1]],1],{i,l-2}];
len=Length[Criticallistmax];
If[deep>len,Print["No hay tantos maximos"]; Abort[];];
tab=Table[Criticallistmax[[i]][[2]]-Min[Take[Criticallistmin,Criticallistmax[[i]][[1]]]],{i,If[deep==0,len,deep]}];
positionmax=Flatten[Position[tab,maxi=Max[{Max[tab],0}]]]//Last;
maxingraph=Criticallistmax[[positionmax]][[1]];
miningraph=Min[minlist=Take[Criticallistmin,maxingraph]];
positionmin=Flatten[Position[minlist,miningraph]]//Last;
Show[ListLinePlot[list,PlotRange->All,PlotStyle->Red,PlotLabel->Style[ToString[maxi]]],ListPlot[{{positionmin+1,list[[positionmin+1]]},{maxingraph+1,list[[maxingraph+1]]}},PlotStyle->{Blue,PointSize[0.015]}],ImageSize->Medium]
];

CharlieMeasureForShowThings[list_]:=CharlieMeasureForShowThings[list,0];

StairCase[x_,eigen_]:=Length[Select[eigen,#<x&]];

StepDecomposition[list_,\[Epsilon]_,elemnts_]:=Select[Flatten[Split[Sort[#],Abs[#1-#2]<\[Epsilon]&]&/@list,1],Length[#]>elemnts&];

(* Rutinas de Cristopher *)

(*Nearest Neighbour Spacings (NNS) of a list. Preferably Apply after Unfold*)
NNS[eigen_]:=Table[Abs[#[[i+1]]-#[[i]]],{i,Length[#]-1}]&[Sort[eigen]];

(*Unfolding of a list*)
Unfold[list_]:=Module[{List0,Staircase,StairTable,FitStaircase,x},List0=Chop[Sort[#]]-Min[#]&[list];Staircase[x_]:=Length[Select[List0,#<x&]];StairTable=Table[{x,Staircase[x]},{x,Min[#],Max[#],Mean[NNS[#]]}]&[List0];FitStaircase = Fit[StairTable, {1,x,x^2,x^3,x^4,x^5,x^6,x^7,x^8,x^9},x];
FitStaircase/.x->#&/@List0];

(*Brody Distribution as function of s, q=0 is Poisson, q=1 is Wigner*)
PBrody[s_,q_]:=Module[{B=Gamma[(2+q)/(1+q)]^(q+1)},B (1+q) s^q*Exp[-B s^(q+1)]];

(*Parametro de Caos LSI, se toma hasta la primera interseccion entre \
Wigner y de Poisson, ya esta region contiene la informacion sobre la \
degeneracion del sistema*)
LSI[unfoldedlist_,bins_] := (4.63551) (Integrate[PDF[HistogramDistribution[unfoldedlist,bins],s],{s,1.0*10^-16,0.472913}] - 0.16109);
(*End of Rutinas de Cristopher*)

(*Amplitude Damping Channel*)
G[t_,\[Lambda]_,\[Omega]0_,\[Gamma]_]:=1/Sqrt[-2 \[Gamma] \[Lambda]+(\[Lambda]+I \[Omega]0)^2] E^(-(1/2) t (\[Lambda]+I \[Omega]0)) (Sqrt[-2 \[Gamma] \[Lambda]+(\[Lambda]+I \[Omega]0)^2] Cosh[1/2 t Sqrt[-2 \[Gamma] \[Lambda]+(\[Lambda]+I \[Omega]0)^2]]+(\[Lambda]+I \[Omega]0) Sinh[1/2 t Sqrt[-2 \[Gamma] \[Lambda]+(\[Lambda]+I \[Omega]0)^2]]);
H[p_]:=-p Log[2,p]-(1-p)Log[2,1-p];
QuantumCapacityDamping[t_]:=If[Abs[G[t,\[Lambda],\[Omega]0,\[Gamma]]]^2>0.5,FindMaximum[H[Abs[G[t,\[Lambda],\[Omega]0,\[Gamma]]]^2 p]-H[(1-Abs[G[t,\[Lambda],\[Omega]0,\[Gamma]]]^2)p],{p,$MinMachineNumber}]//First,0];
ClassicalCapacityDamping[t_]:=FindMaximum[H[p]+H[Abs[G[t,\[Lambda],\[Omega]0,\[Gamma]]]^2 p]-H[(1-Abs[G[t,\[Lambda],\[Omega]0,\[Gamma]]]^2)p],{p,$MinMachineNumber},MaxIterations->Infinity]//First;
(*AveK[list_]:=(Last[list][[1]]-First[list][[1]])^(-1)(list[[2]][[1]]-list[[1]][[1]])Sum[list[[All,2]][[i]],{i,Length[list]}];*)

(*Routines for check divisibility properties*)
BasisElement[i_,j_]:=Table[If[k==i&&l==j,1,0],{k,2},{l,2}];
BasisElementOneIndex[i_]:=Switch[i,1,BasisElement[1,1],2,BasisElement[1,2],3,BasisElement[2,1],4,BasisElement[2,2]];
w=Table[Tr[BasisElementOneIndex[i+1].PauliMatrix[j]/Sqrt[2]],{i,0,3},{j,0,3}];\[Omega]=Proyector[Bell[2]];\[Omega]ort=IdentityMatrix[4]-\[Omega];

DivisivilityKindOf[\[Lambda]1_,\[Lambda]2_,\[Lambda]3_]:=Module[{eigen,list,L},
list=Sqrt[Sort[{1,\[Lambda]1^2,\[Lambda]2^2,\[Lambda]3^2}]];
L=MatrixLog[Chop[w.{{1,0,0,0},{0,\[Lambda]1,0,0},{0,0,\[Lambda]2,0},{0,0,0,\[Lambda]3}}.Dagger[w]]];
If[
(*Checking Complete Positivity*)
1+\[Lambda]1+\[Lambda]2+\[Lambda]3>=0&&1-\[Lambda]1-\[Lambda]2+\[Lambda]3>=0&&1-\[Lambda]1+\[Lambda]2-\[Lambda]3>=0&&1+\[Lambda]1-\[Lambda]2-\[Lambda]3>=0,
If[
(*Evaluating CP-divisibility and p-divisibility*)
(*Evaluating for p-divsibility*)\[Lambda]1 \[Lambda]2 \[Lambda]3>0,
If[ (*Evaluating for CP-div*)
list[[1]]^2*list[[4]]^2>=Product[list[[i]],{i,1,4}],
(*Evaluating for markov type evolution*)
If[
PositiveSemidefiniteMatrixQ[\[Omega]ort.Reshuffle[L].\[Omega]ort]&&HermitianMatrixQ[L],4,3
],2],1],0]]
End[] 
EndPackage[]
