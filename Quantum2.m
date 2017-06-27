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
CharlieMeasureAveForShowThings::usage = "CharlieMeasureAveForShowThings[list_,deep_]"
StairCase::usage = "StairCase[x_,eigen_]"
NNS::usage = "NNS[eigen_]"
Unfold::usage = "Unfold[list_]"
PBrody::usage = "PBrody[s_,q_]"
LSI::usage = "LSI[unfoldedlist_,bins_] Level Statistic Indicator or simply Gamma parameter"
G::usage = "Amplitude Damping Quantum Channel G[t_,\[Lambda]_,\[Omega]0_,\[Gamma]_]"
H::usage = "Binary Shannon Entropy H[p_]"
QuantumCapacityDamping::usage = "Quantum Capacity of the quantum damping, It must be specified the parameters of G, QuantumCapacityDamping[t_,\[Lambda]_,\[Gamma]_,\[Omega]0_]"
ClassicalCapacityDamping::usage = "EA Classical Capacity of the quantum damping, It must be specified the parameters of G, ClassicalCapacityDamping[t_,\[Lambda]_,\[Gamma]_,\[Omega]0_]"
StepDecomposition::usage = "StepDecomposition[list_,\[Epsilon]_,elemnts_]"
BasisElement::usage = "BasisElement[i_,j_] Dont worry about this, not yet"
BasisElementOneIndex::usage = "BasisElementOneIndex[i_] Dont worry about this, not yet"
DivisibilityKindOf::usage = "DivisivilityKindOf[\[Lambda]1_,\[Lambda]2_,\[Lambda]3_] The lambdas state for the singular values of the unital channel
up rotations), then this function gives 0 when the channel is not CPTP, 1 if the channel is CPTP, 2 if it is p-divisible, 3 if it is compatible
with CP-divisible dynamics and 4 if the channel can be written as exp(L) with L Lindblad."
EntanglementBreakingQ::usage = "EntanglementBreakingQ[x_,y_,z_] this function checks if the channel is entanglement-breaking 
in the sense of a separable Jamilokowski state."
DivisibilityKindOfGeneral::usage = "DivisibilityKindOfGeneral[channel_], Sure it works at least for channels with diagonal lorentz form in the case of CP-divisibility"
gRHP::usage = "gRHP[list_] Calculation of the Rivas g(t) from fidelity, i. e. from D(t) for dephasing channels."
PositiveDerivatives::usage = "PositiveDerivatives[list_] etc."
maximizer::usage = "maximizer[list_] divides the second column of the list by the maximum value of the original list."
QuantumMapInPauliBasis::usage = "QuantumMapInPauliBasis[channel_] This function constructs the Pauli basis channel representation of one qubit"
QuantumMapInUnitBasis::usage = "QuantumMapInInUnitBasis[channel_] This function constructs the Pauli basis channel representation of one qubit"
FromPauliToUnit::usage = " Takes channels in Pauli basis to unital matrices basis."
FromUnitToPauli::usage = " Oposite of From Unit to Pauli."
UnitalChannelInPauliBasis::usage = "UnitalChannelInPauliBasis[x_,y_,z_], where x,y and z are the weights of the corresponding unitaries, by convexity teh weight of the identity operation is automatically determined "
EigensystemOrdered::usage = "EigensystemOrdered[matrix_], Routine for compute Eigensystem with ordered Eigenvalues, see code for further information, this is useful for several routines which require the same basis ordering for different diagonalizations."
RealMatrixLogarithmComplexCase::usage = "RealMatrixLogarithmComplexCase[matrix,k=0] Real logarithm of a matrix for the complex eigenvalues case, returns Falese if the logarithm doesnt exists."
RealMatrixLogarithmRealCase::usage = "RealMatrixLogarithmComplexCase[matrix,k=0] Real logarithm of a matrix for the real eigenvalues case, returns Falese if the logarithm doesnt exists."
HermiticityPreservingAndCCPOfGenerator::usage= "HermiticityPreservingAndCCPOfGenerator[matrix,upto_Branch] Test if the given channels has a hermiticity preserving and ccp generator, returns False if there is no such generator."
HasHermitianPreservingAndCCPGenerator::usage = "HasHermiticitianPreservingAndCCPOfGenerator[matrix_,upto_Branch] Returns True if the channels has hermitian preserving and ccp generator."
DecompositionOfChannelsInSO31::usage = "Performs decomposition in orthochronus Lorentz group of the matrix, returns False if the decomposition can not be done."
LorentzMatrixQ::usage = "LorentzMatrixQ[matrix_] returns if the given matrix is a valid Lorentz transformation with the signature +---."
ForceSameSignatureForLorentz::usage = "ForceSameSignatureForLorentz[matrix_] Forces or fixes the signature of the given Lorentz transformation."
QuantumMaptoR::usage = "QuantumMaptoR[channel_] Representation of Jamiolokowski state in \[Sigma]iotimes \[Sigma]j basis."
DiagonalMatrixQ::usage = "DiagonalMatrixQ[matrix_] Works similar to the other matrix tests of mathematica."
PositiveSemidefiniteMatrixCustomQ::usage = "Mathematica test gives strange results, this routine check such test hunred times."
PositiveSemidefiniteMatrixCustom2Q::usage = "Custom check."
PositiveSemidefiniteMatrixCustom3Q::usage = "By eigenvalues of the hermitian part"
DecompositionOfChannelsInSO::usage = "Singular value decomposition but in SO."
HermitianPart::usage = "HermitianPart."
DecompositionOfChannelsInSO31Diagonal::usage= "The same of DecompositionOfChannelsInSO31 but forcing a diagonal in the spatial block."
AddOne::usage="Direct sum from behind with an identity of dimension one."
InfinitySupressor::usage = "InfinitySupressor[lista_,infinity_] This functions makes infinities finide by a given bound."
PartialDecompositionofChannelInSO::usage = "PartialDecompositionofChannelInSO[matrix_]."
TestViaRankOfCPDIV::usage = "TestViaRankOfCPDIV[matrix_], gives True or False."
KrausRank::usage = "Computes the Kraus rank of qubit channels given in Pauli basis."
WolfEisertCubittCiracMeasureQubitCase::usage = "WolfEisertCubittCiracMeasureQubitCase[channel_] Computes waht the name says, checking up to 10 branches."
HilbertSpaceBasisParametrization::usage = "HilbertSpaceBasisParametrization[n_] A most general parametrization of the basis of C^n, or equivalently an element of SU(n)."
PhasesEliminator::usage = "PhasesEliminator[basis_] Removes global phases of the basis given by HilbertSpaceBasisParametrization."


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

CharlieMeasure[list_] := 
  Module[{l, listD, Criticallistmin, Criticallistmax, position, maxi, 
    len},
   l = Length[list];
   listD = Table[list[[i + 1]] - list[[i]], {i, l - 1}];
   Criticallistmax = 
    Sort[DeleteCases[
      Table[If[(listD[[i]] > 0 && 
           listD[[i + 1]] <= 0) || (listD[[i]] >= 0 && 
           listD[[i + 1]] < 0), {i, list[[i + 1]]}, 0], {i, l - 2}], 
      0], #1[[2]] > #2[[2]] &];
len=Length[Criticallistmax];
If[len==0,0,Max[{0,Max[Table[Criticallistmax[[i]][[2]]-Min[Take[list,Criticallistmax[[i]][[1]]]],{i,1,len}]]}]]
];

CharlieMeasureAve[list_] := 
  Module[{l, listD, Criticallistmin, Criticallistmax, position, maxi, 
    len},
   l = Length[list];
   listD = Table[list[[i + 1]] - list[[i]], {i, l - 1}];
   Criticallistmax = 
    Sort[DeleteCases[
      Table[If[(listD[[i]] > 0 && 
           listD[[i + 1]] <= 0) || (listD[[i]] >= 0 && 
           listD[[i + 1]] < 0), {i, list[[i + 1]]}, 0], {i, l - 2}], 
      0], #1[[2]] > #2[[2]] &];
len=Length[Criticallistmax];
If[len==0,0,Max[{0,Max[Table[Criticallistmax[[i]][[2]]-Mean[Take[list,Criticallistmax[[i]][[1]]-1]],{i,1,len}]]}]]
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

CharlieMeasureAveForShowThings[list_,deep_]:=Module[{mean,l,listD,Criticallistmin,Criticallistmax,position,maxi,len,tab,positionmax,miningraph,maxingraph,minlist,positionmin},
l=Length[list];
listD=Table[list[[i+1]]-list[[i]],{i,l-1}];
Criticallistmax=Sort[DeleteCases[Table[If[(listD[[i]]>0&&listD[[i+1]]<=0)||(listD[[i]]>=0&&listD[[i+1]]<0),{i,list[[i+1]]},0],{i,l-2}],0],#1[[2]]>#2[[2]]&];
Criticallistmin=Table[If[(listD[[i]]<=0&&listD[[i+1]]>0)||(listD[[i]]<0&&listD[[i+1]]>=0),list[[i+1]],1],{i,l-2}];
len=Length[Criticallistmax];
If[deep>len,Print["No hay tantos maximos"]; Abort[];];
tab=Table[Criticallistmax[[i]][[2]]-Mean[Take[list,Criticallistmax[[i]][[1]]-1]],{i,If[deep==0,len,deep]}];
positionmax=Flatten[Position[tab,maxi=Max[{Max[tab],0}]]]//Last;
maxingraph=Criticallistmax[[positionmax]][[1]];
mean=Mean[Take[list,Criticallistmax[[positionmax]][[1]]-1]];
Show[ListLinePlot[list,PlotRange->All,PlotStyle->Red,PlotLabel->Style[ToString[maxi]]],ListPlot[{{maxingraph+1,list[[maxingraph+1]]}},PlotStyle->{Blue,PointSize[0.015]}],Plot[mean,{x,0,maxingraph},PlotStyle->Directive[Black,Thick]],ImageSize->Medium]
];

CharlieMeasureAveForShowThings[list_]:=CharlieMeasureAveForShowThings[list,0];

CharlieMeasureForShowThings[list_]:=CharlieMeasureForShowThings[list,0];

StairCase[x_,eigen_]:=Length[Select[eigen,#<x&]];

StepDecomposition[list_,\[Epsilon]_,elemnts_]:=Select[Flatten[Split[Sort[#],Abs[#1-#2]<\[Epsilon]&]&/@list,1],Length[#]>elemnts&];

(* Rutinas de Cristopher *)

(*Nearest Neighbour Spacings (NNS) of a list. Preferably Apply after Unfold*)
NNS[eigen_]:=Table[Abs[#[[i+1]]-#[[i]]],{i,Length[#]-1}]&[Sort[eigen]];

(*Unfolding of a list*)
Unfold[list_]:=Module[{List0,Staircase,StairTable,FitStaircase,x},List0=Chop[Sort[#]]-Min[#]&[list];Staircase[x_]:=Length[Select[List0,#<x&]];StairTable=Table[{x,Staircase[x]},{x,Min[#],Max[#],Mean[NNS[#]]}]&[List0];FitStaircase = Fit[StairTable, {1,x,x^2,x^3,x^4,x^5,x^6,x^7,x^8,x^9},x];
FitStaircase/.x->#&/@List0];

Unfold[list_, grado_] :=
  Block[{List0, Staircase, StairTable, FitStaircase, x},
   List0 = Chop[Sort[#]] - Min[#] &[list];
   Staircase[x_] := Length[Select[List0, # < x &]];
   StairTable =
    Table[{x, Staircase[x]}, {x, Min[#], Max[#], Mean[NNS[#]]}] &[
     List0]; FitStaircase =
    Fit[StairTable, Table[x^n, {n, 0, grado}], x];
   FitStaircase /. x -> # & /@ List0];

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
QuantumCapacityDamping[t_,\[Lambda]_,\[Gamma]_,\[Omega]0_]:=If[Abs[G[t,\[Lambda],\[Omega]0,\[Gamma]]]^2>0.5,FindMaximum[H[Abs[G[t,\[Lambda],\[Omega]0,\[Gamma]]]^2 p]-H[(1-Abs[G[t,\[Lambda],\[Omega]0,\[Gamma]]]^2)p],{p,$MinMachineNumber}]//First,0];
ClassicalCapacityDamping[t_,\[Lambda]_,\[Gamma]_,\[Omega]0_]:=FindMaximum[H[p]+H[Abs[G[t,\[Lambda],\[Omega]0,\[Gamma]]]^2 p]-H[(1-Abs[G[t,\[Lambda],\[Omega]0,\[Gamma]]]^2)p],{p,$MinMachineNumber},MaxIterations->Infinity]//First;
(*AveK[list_]:=(Last[list][[1]]-First[list][[1]])^(-1)(list[[2]][[1]]-list[[1]][[1]])Sum[list[[All,2]][[i]],{i,Length[list]}];*)

(*Routines for check divisibility properties*)
BasisElement[i_,j_]:=Table[If[k==i&&l==j,1,0],{k,2},{l,2}];
BasisElementOneIndex[i_]:=Switch[i,1,BasisElement[1,1],2,BasisElement[1,2],3,BasisElement[2,1],4,BasisElement[2,2]];
w=Table[Tr[Dagger[BasisElementOneIndex[i+1]].PauliMatrix[j]/Sqrt[2]],{i,0,3},{j,0,3}];\[Omega]=Proyector[Bell[2]];\[Omega]ort=IdentityMatrix[4]-\[Omega];

DivisibilityKindOf[\[Lambda]1_,\[Lambda]2_,\[Lambda]3_]:=Module[{eigen,list},
list=Sort[Sqrt[{1,\[Lambda]1^2,\[Lambda]2^2,\[Lambda]3^2}]];
If[
(*Checking Complete Positivity*)
Chop[1+\[Lambda]1+\[Lambda]2+\[Lambda]3]>=0&&Chop[1-\[Lambda]1-\[Lambda]2+\[Lambda]3]>=0&&Chop[1-\[Lambda]1+\[Lambda]2-\[Lambda]3]>=0&&Chop[1+\[Lambda]1-\[Lambda]2-\[Lambda]3]>=0,
If[
(*Evaluating CP-divisibility and p-divisibility*)
(*Evaluating for p-divsibility*)\[Lambda]1 \[Lambda]2 \[Lambda]3>=0,
If[ (*Evaluating for CP-div*)
list[[1]]^2>=\[Lambda]1*\[Lambda]2*\[Lambda]3,
(*Evaluating for markov type evolution*)
If[(*checking hermiticity preserving and CCP*)
HermiticityPreservingAndCCPOfTheGeneratorQForDiagonal[\[Lambda]1,\[Lambda]2,\[Lambda]3,2]
,4,3
],2],1],0]];

DivisibilityKindOf[\[Lambda]_]:=DivisibilityKindOf[\[Lambda][[1]],\[Lambda][[2]],\[Lambda][[3]]];

PositiveSemidefiniteMatrixCustomQ[matrix_]:=And@@Table[PositiveSemidefiniteMatrixQ[matrix],{1000}];

PositiveSemidefiniteMatrixCustom2Q[matrix_]:=And@@Table[v=RandomState[Length[matrix]];Chop[QuantumDotProduct[v,HermitianPart[matrix].v]]>=0,{1000}];

PositiveSemidefiniteMatrixCustom3Q[matrix_]:=And@@NonNegative[Chop[Eigenvalues[HermitianPart[matrix]]]];

g=DiagonalMatrix[{1,-1,-1,-1}];

TestViaRankOfCPDIV[matrix_]:=Module[{c,rank,var,det,list},
c=matrix.g.Transpose[matrix].g//Chop;
rank=MatrixRank[matrix];
det=Det[matrix];
list=Sqrt[Abs[Eigenvalues[c]]];
If[((var=DiagonalizableMatrixQ[matrix])||(MatrixRank[c]<rank))&&det>=0,
If[var&&(rank<2||(Chop[Apply[Times,Sort[list^2][[{1,4}]]]]>=(Times@@list)&&det>0)),True,
If[(MatrixRank[c]<rank),True,False
],False]
,False]
];


DivisibilityKindOfGeneral[channel_,branches_:1]:=Module[{eigen,list,tmp,det,form},
If[
(*Checking Complete Positivity*)
PositiveSemidefiniteMatrixCustom3Q[tmp=Reshuffle[Chop[FromPauliToUnit[Chop[channel]]]]],
If[
HasHermitianPreservingAndCCPGenerator[Chop[channel],branches],4,
If[TestViaRankOfCPDIV[channel],3,
If[Chop[Det[channel]]>0,2,1]]
],0]*0.25
];

EntanglementBreakingQ[x_,y_,z_]:=If[DivisibilityKindOf[x,y,z]>0,If[Max[0,1/4 (-Abs[-1+x+y-z]-Abs[-1+x-y+z]-Abs[-1-x+y+z]-Abs[1+x+y+z]+8 Max[1/4 Abs[-1+x+y-z],1/4 Abs[-1+x-y+z],1/4 Abs[-1-x+y+z],1/4 Abs[1+x+y+z]])]<=0,2,1],0];

EntanglementBreakingQ[channel_]:=If[Concurrence[0.5*Reshuffle[FromPauliToUnit[channel]]]==0,True,False];

gRHP[list_]:=Map[If[#<0,0,#]&,Table[list[[i+1]]/list[[i]]-1,{i,Length[list]-1}]];
PositiveDerivatives[list_]:=Map[If[#<0,0,#]&,Table[list[[i+1]]-list[[i]],{i,Length[list]-1}]];
(*Misc*)

maximizer[list_,factor_]:=Module[{max},
max=Max[list[[All,2]]];
Map[{#[[1]],#[[2]]/(factor*max)}&,list]
];

maximizer[list_,factor_,deep_]:=Module[{max,maxlist,pos,listnew},
maxlist=Sort[list[[All,2]],Greater];
pos=Table[Position[list[[All,2]],maxlist[[i]]],{i,deep}]//Transpose//First;
listnew=Delete[list,pos];
max=Max[listnew[[All,2]]];
Map[{#[[1]],#[[2]]/(factor*max)}&,listnew]
]

maximizer[list_,factor_,0]:=maximizer[list,factor];

maximizer[list_]:=maximizer[list,1,0];

QuantumMapInPauliBasis[channel_]:=1/2Table[Tr[PauliMatrix[i].channel[PauliMatrix[j]]],{i,0,3},{j,0,3}];

QuantumMapInUnitBasis[channel_]:=Table[Tr[Dagger[BasisElementOneIndex[i]].channel[BasisElementOneIndex[j]]],{i,1,4},{j,1,4}];

FromPauliToUnit[channel_]:=w.channel.Dagger[w];
FromUnitToPauli[channel_]:=Dagger[w].channel.w;

EigensystemOrdered[A_]:=Module[{eigs,vecs,list1,list2},
{eigs,vecs}=Eigensystem[A];
list1=Partition[Riffle[eigs,vecs],2];
list2=Sort[Sort[list1,Re[#1[[1]]]<Re[#2[[1]]]&],Im[#1[[1]]]>Im[#2[[1]]]&];
list2//Transpose
];

RealMatrixLogarithmRealCase[matrix_,k_:0]:=Module[{eig,eigneg,vectors,V,L,pos},
{eig,vectors}=Eigensystem[matrix];
If[Element[eig,Reals]==False,Print["Se uso la rutina para logaritmo complejo"];RealMatrixLogarithmComplexCase[matrix,k],
V=Transpose[vectors]//Chop;
eigneg=Select[eig,#<0&]//Chop;
L=DiagonalMatrix[Log[Abs[eig]]]//Chop;
If[Length[eigneg]==0,V.L.Inverse[V],
If[Length[eigneg]==2&&Chop[eigneg[[1]]-eigneg[[2]],10^(-10)]==0,
pos=Position[eig,eigneg[[1]]][[{1,2},1]];
L[[pos[[1]],pos[[2]]]]=(2 k+1) Pi;
L[[pos[[2]],pos[[1]]]]=-(2 k+1) Pi;
V.L.Inverse[V]//Chop
,False]]
]
];

RealMatrixLogarithmComplexCase[channel_,k_:0]:=Module[{mat,diag,w,e,pos,list,d2,w2,chreorlog,wreal,list2,is},
{diag,w}=Chop[EigensystemOrdered[channel]]//Chop;
If[Element[diag,Reals],(*Print["Se uso la rutina para logaritmo real"];*)RealMatrixLogarithmRealCase[channel,k],
mat=DiagonalMatrix[diag];
w=Transpose[w];
e=ConstantArray[0.0,Dimensions[channel]];
list=Select[diag,Chop[Im[#]]!=0&]//Chop;
list2=Select[diag,Chop[Im[#]]==0&]//Chop;
If[(Length[list2]==2&&AllTrue[list2,NonNegative])==False||(Length[list2]==1)==True,is=False,is=True];
pos=Flatten[Table[Position[diag,i],{i,list}]];
mat[[pos[[1]],pos[[1]]]]=Re[diag[[pos[[1]]]]];
mat[[pos[[2]],pos[[2]]]]=Re[diag[[pos[[1]]]]];
mat[[pos[[1]],pos[[2]]]]=Im[diag[[pos[[1]]]]];
mat[[pos[[2]],pos[[1]]]]=-Im[diag[[pos[[1]]]]];
{d2,w2}=EigensystemOrdered[mat]//Chop;
e[[pos[[1]],pos[[2]]]]=1;
e[[pos[[2]],pos[[1]]]]=-1;
chreorlog=Chop[MatrixLog[mat]+2 Pi k e];
If[Total[Flatten[Chop[d2-diag]]]==0,

w2=Transpose[w2]//Chop;
d2=DiagonalMatrix[d2]//Chop;
wreal=w.Inverse[w2]//Chop;
If[is,wreal.chreorlog.Inverse[wreal]//Chop,is]

,Return["bad calculation"]
]
]
];


HermiticityPreservingAndCCPOfGenerator[matrix_,upto_:1]:=Module[{is,L,i,branches,b},
b=0;
branches=Table[b=b+(-1)^(j+1)*j,{j,0,2*upto}];
is=False;
If[DiagonalizableMatrixQ[matrix],
Table[If[PositiveSemidefiniteMatrixCustom3Q[Chop[\[Omega]ort.Chop[Reshuffle[FromPauliToUnit[L=RealMatrixLogarithmComplexCase[Chop[matrix],k]//Chop]]].\[Omega]ort]],is=True;i=k;Return[Null,Table],is=False;],{k,branches}];,
Return["non diagonalizable"];
];
If[i!=0,Print["El logaritmo es real hasta k= "<>ToString[i]]];
If[is==True,L,False]
];

DecompositionOfChannelsInSO[map_]:=Module[{a,e,i,newmap,eig,vecs},
{a,e,i}=SingularValueDecomposition[map];
{eig,vecs}=Chop[Eigensystem[Det[a]*Det[i]*(i.e.Transpose[i])]];
{a.Transpose[i].Transpose[vecs],DiagonalMatrix[eig],vecs}
];

ForceSameSignatureForLorentz[matrix_]:=Module[{tmpmat,eigs,o,tildeM},
tmpmat=matrix.g.Transpose[matrix]//Chop;
{eigs,o}=EigensystemOrdered[tmpmat];tildeM=DiagonalMatrix[eigs];
{o.matrix,o}
];

HermitianPart[m_]:=0.5*(m+Dagger[m]);

QuantumMaptoR[channel_]:=Module[{g,\[Rho]},
g=DiagonalMatrix[{1,-1,-1,-1}];
\[Rho]=Reshuffle[FromPauliToUnit[channel]]/2//Chop;
Table[Tr[\[Rho].KroneckerProduct[PauliMatrix[i],PauliMatrix[j]]],{i,0,3},{j,0,3}]//FullSimplify//Chop
];

LorentzMatrixQ[matrix_]:=
Chop[matrix[[1,1]]]>0&&Chop[Det[matrix]]==1&&Chop[matrix.g.Transpose[matrix]]==g;

DecompositionOfChannelsInSO31[matrix_]:=Module[{c,x,j,n,eig1,eig2,o,d,leftL,rightL},
c=g.matrix.g.Transpose[matrix]//Chop;
{x,j}=SchurDecomposition[c]//Chop;x=Inverse[x]//Chop;
n=x.g.Transpose[x]//Chop;
{eig1,o}=EigensystemOrdered[n]//Chop;
leftL=Transpose[x].Transpose[o]//Chop;
c=g.Transpose[matrix].g.matrix;
{x,j}=SchurDecomposition[c]//Chop;x=Inverse[x]//Chop;
n=x.g.Transpose[x]//Chop;
{eig2,o}=EigensystemOrdered[n]//Chop;
If[Chop[eig1]==Diagonal[g]&&eig2==Diagonal[g],None,Print["Wrong calculation"];Abort[];];
rightL=Transpose[x].Transpose[o]//Chop;
If[Det[rightL]==-1,rightL=Sign[rightL[[1,1]]]g.rightL;,If[rightL[[1,1]]<0,rightL=-g.g.rightL;]];
If[Det[leftL]==-1,leftL=Sign[leftL[[1,1]]]g.leftL;,If[leftL[[1,1]]<0,leftL=-g.g.leftL;]];
If[(LorentzMatrixQ[leftL]&&LorentzMatrixQ[rightL])==False,Print["Decomposition not done"];,{leftL,Transpose[leftL].matrix.rightL,rightL}//Chop]
];

AddOne[matrix_]:=Flatten[{{{1.0,0.0,0.0,0.0}},Table[Flatten[{0.0,matrix[[i,All]]}],{i,1,3}]},1];

DecompositionOfChannelsInSO31Diagonal[matrix_]:=Module[{c,x,j,n,eig1,eig2,o,d,leftL,rightL,form,aux,form2,a,e,i},
(*Transformada Izquierda:*)
c=g.matrix.g.Transpose[matrix]//Chop;
{x,j}=SchurDecomposition[c]//Chop;x=Inverse[x]//Chop;
n=x.g.Transpose[x]//Chop;
{eig1,o}=EigensystemOrdered[n]//Chop;
leftL=Transpose[x].Transpose[o]//Chop;
(*Transformada Derecha:*)
c=g.Transpose[matrix].g.matrix;
{x,j}=SchurDecomposition[c]//Chop;x=Inverse[x]//Chop;
n=x.g.Transpose[x]//Chop;
{eig2,o}=EigensystemOrdered[n]//Chop;
If[Chop[eig1]==Diagonal[g]&&eig2==Diagonal[g],None,Print["Wrong calculation"];Abort[];];
rightL=Transpose[x].Transpose[o]//Chop;
(*Se define la forma normal preliminar, pues puede que no salga diagonal*)
form=Transpose[leftL].matrix.rightL;
e=form;
(*En este paso construyo un canal solo con el bloque espacial de la matriz: *)
form2=AddOne[Take[form,{2,4},{2,4}]];
(*Por si no sale diagonal el bloque espacial, aqui se fuerza descomponiendolo en SO(3):*)
If[DiagonalMatrixQ[form2]==False,
{a,e,i}=DecompositionOfChannelsInSO[form2]//Chop;
aux=form-form2;
d=aux[[All,1]];
d=Transpose[a].d//Chop;
e[[All,1]]=d;e[[1,1]]=1;
(*Se redefinen las transformaciones de Lorentz con las nuevas rotaciones:*)
leftL=leftL.a//Chop;
rightL=rightL.Transpose[i]//Chop;
];
(*Se hacen propias y ortocronas usando el grupo cociente O(3,1)/SO^+(3,1):*)
If[Det[rightL]==-1,rightL=Sign[rightL[[1,1]]]g.rightL;,If[rightL[[1,1]]<0,rightL=-g.g.rightL;]];
If[Det[leftL]==-1,leftL=Sign[leftL[[1,1]]]g.leftL;,If[leftL[[1,1]]<0,leftL=-g.g.leftL;]];
(*Un test util para mantener esta funcion y el output:*)
If[(LorentzMatrixQ[leftL]&&LorentzMatrixQ[rightL])==False,Print["Decomposition not done"];,{leftL,e,rightL}//Chop]
];

PartialDecompositionofChannelInSO[matrix_]:=Module[{form2,a,e,i,aux,form,d,leftL,rightL},
form2=AddOne[Take[matrix,{2,4},{2,4}]];
{a,e,i}=DecompositionOfChannelsInSO[form2]//Chop;
aux=matrix-form2;
d=aux[[All,1]];
d=Transpose[a].d//Chop;
e[[All,1]]=d;e[[1,1]]=1;
leftL=a//Chop;
rightL=Transpose[i]//Chop;
{leftL,e,rightL}
];

HasHermitianPreservingAndCCPGenerator[matrix_,upto_:1]:=Module[{is,L,i,branches,b},
b=0;
branches=Table[b=b+(-1)^(j+1)*j,{j,0,2*upto}];
is=False;
If[DiagonalizableMatrixQ[matrix],

Table[If[PositiveSemidefiniteMatrixCustom3Q[\[Omega]ort.FullSimplify[Reshuffle[FromPauliToUnit[L=RealMatrixLogarithmComplexCase[Chop[matrix],k]//Chop]//Chop]].\[Omega]ort//Chop],is=True;i=k;Return[Null,Table],is=False;],{k,branches}];,
Return["non diagonalizable"];
];
If[i!=0,Print["Hermiticity preserving and ccp condition is fulfilled until k= "<>ToString[i]]];
If[is==True,is,False]
];

WolfEisertCubittCiracMeasureQubitCase[channel_,OptionsPattern[{branches->5,noisedelta->0.01}]]:=Module[{dev,noise,lol,\[Mu]},
noise[\[Mu]_]:=-\[Mu] DiagonalMatrix[{0,1,1,1}];
dev=RealMatrixLogarithmComplexCase[channel,10];
If[Length[dev]==0,0,
lol=Min[DeleteCases[Table[If[PositiveSemidefiniteMatrixCustom3Q[Chop[\[Omega]ort.Reshuffle[FromPauliToUnit[RealMatrixLogarithmComplexCase[channel,i]+noise[j]]].\[Omega]ort]],j,None],{i,0,OptionValue[branches]},{j,0.0,1.0,OptionValue[noisedelta]}]//Flatten,None]];
Exp[-3lol]
]
];

UnitalChannelInPauliBasis[x_,y_,z_]:=DiagonalMatrix[(1-x-y-z){1,1,1,1}+x{1,1,-1,-1}+y{1,-1,1,-1}+z{1,-1,-1,1}];
UnitalChannelInPauliBasis[vec_]:=UnitalChannelInPauliBasis[vec[[1]],vec[[2]],vec[[3]]];

DiagonalMatrixQ[matrix_]:=If[Total[Abs[Flatten[DiagonalMatrix[Diagonal[matrix]]-matrix]]]==0,True,False];

InfinitySupressor[lista_,infinity_]:=Module[{list},
list=lista;
Table[If[list[[i]][[2]]>infinity,list[[i]][[2]]=infinity],{i,Length[list]}];
list
];

KrausRank[channel_]:=MatrixRank[Reshuffle[FromPauliToUnit[channel]]];
	
End[]

SinCosList[j_,vec_]:=Table[If[i==j,Cos,Sin][\[Theta][vec][i+1]],{i,0,j}];	
	
SinList[j_,vec_]:=Table[Sin[\[Theta][vec][i+1]],{i,0,j}];

SinTerm[j_,vec_]:=Times@@SinList[j,vec];

SinCosTerm[j_,vec_]:=Times@@SinCosList[j,vec];

Vector[index_,dim_,basis_]:=Sum[If[dim-1==i,SinTerm[i-1,index],SinCosTerm[i,index]]Exp[I \[Phi][index][i]]basis[[i+1]],{i,0,dim-1}];

Vector[index_,dim_]:=Table[If[dim-1==i,SinTerm[i-1,index],SinCosTerm[i,index]]Exp[I \[Phi][index][i]],{i,0,dim-1}];

vecuptoU[vec_,index_,b_]:=ReplaceAll[D[vec,\[Theta][index][b]],Table[\[Theta][index][i]->Pi/2,{i,1,b-1}]];

BasisComplement[vec_,index_,dim_]:=Table[vecuptoU[vec,index,i],{i,dim-(index+1)}];

HilbertSpaceBasisParametrization[n_]:=Module[{vec,basis},
vec=Vector[0,n];
Flatten[{{vec},Table[
basis=BasisComplement[vec,i-1,n];
vec=Vector[i,n-i,basis];
vec
,{i,n-1}]},1]
];

PhasesEliminator[basis_]:=Module[{len},
len=Length[basis];
ReplaceAll[basis,Table[\[Phi][i][0]->0,{i,0,len-1}]]
];
 
EndPackage[]



