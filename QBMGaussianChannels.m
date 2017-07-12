(* ::Package:: *)

BeginPackage["QBMGaussianChannels`",{"Quantum2`"}]
matrixTstandAlone::usage = "matrixTstandAlone[A_,DA_,DDA_] or matrixTstandAlone[t_]";
matrixNstandAlone::usage = "matrixNstandAlone[A_,DA_,DDA_,S_,DS_,DDS_,qsq_,psq_] or matrixNstandAlone[t_]";
qsq::usage = "Just call it." 
psq::usage = "Just call it." 
GLaplace::usage = "GLaplace[z_]."
GFourier::usage = "GFourier[\[Omega]_]."
SDesdeFourierNumerico::usage = "SDesdeFourierNumerico[t_]."
SDesdeFourierAnalitico::usage = "SDesdeFourierAnalitico[\[Tau]_]."
ADesdeLaplace::usage = "ADesdeLaplace[t_]."
ADesdeFourier::usage = "ADesdeFourier[t_]."
A::usage = "A[t_]."
DA::usage = "DA[t_]."
DDA::usage = "DDA[t_]."
S::usage = "S[t_]."
DS::usage = "DS[t_]."
DDS::usage = "DDS[t_]."
\[Gamma]Drude::usage = "\[Gamma]Drude[\[Omega]_] Damping Kernel for Drude in Fourier Representation."
CompositionLaw::usage = "CompositionLaw[A_,B_]."
InverseTN::usage = "InverseTN[A_]."
cptpCondition::usage = "cptpCondition[A_]."
cptpMatrix::usage = "cptpMatrix[A_]."
ComputeF::usage = "ComputeF[A_]. Computes the Illuminatti indicator of divisibility."
ComputeFlist::usage = "ComputeFlist[listcorr_]."
ComputecptpMatrixlist::usage = "ComputecptpMatrixlist[listcorr_]."
ComputecptpMatrixComponents::usage = "ComputecptpMatrixComponents[listcorr_,which_]."


(* ::Input::Initialization:: *)
Begin["Private`"] 
\[CapitalOmega]={{0,1},{-1,0}};
matrixTstandAlone[A_,DA_,DDA_]:={{-((2 M DA)/hbar),(2 A)/hbar},{(2 M^2 DDA)/hbar,-((2 M DA)/hbar)}}//Chop;
matrixNstandAlone[A_,DA_,DDA_,S_,DS_,DDS_,qsq_,psq_]:={{qsq+(4 A (psq A-hbar M DS))/hbar^2-S^2/qsq,-((2 M (A (-2 psq DA+hbar M DDS)+hbar M DA DS))/hbar^2)-(M DS S)/qsq},{-((2 M (A (-2 psq DA+hbar M DDS)+hbar M DA DS))/hbar^2)-(M DS S)/qsq,psq+M^2 ((4 DA (psq DA-hbar M DDS))/hbar^2-DS^2/qsq)}}//Chop;
matrixTstandAlone[t_]:={{-((2 M DA[t])/hbar),(2 A[t])/hbar},{(2 M^2 DDA[t])/hbar,-((2 M DA[t])/hbar)}};
matrixNstandAlone[t_]:={{qsq+(4 A[t] (psq A[t]-hbar M DS[t]))/hbar^2-S[t]^2/qsq,-((2 M (A[t] (-2 psq DA[t]+hbar M DDS[t])+hbar M DA[t] DS[t]))/hbar^2)-(M DS[t] S[t])/qsq},{-((2 M (A[t] (-2 psq DA[t]+hbar M DDS[t])+hbar M DA[t] DS[t]))/hbar^2)-(M DS[t] S[t])/qsq,psq+M^2 ((4 DA[t] (psq DA[t]-hbar M DDS[t]))/hbar^2-DS[t]^2/qsq)}};
\[HBar]:=hbar;

(*Desviaciones*)qsq:=1/(M \[Beta]\[Beta]) Sum[1/(\[Omega]0^2+(2Pi n/(\[Beta]\[Beta] hbar))^2+Abs[2Pi n/(\[Beta]\[Beta] hbar)] \[Gamma]hat[Abs[2Pi n/(\[Beta]\[Beta] hbar)]]),{n,-1000,1000}];
psq:=M/\[Beta]\[Beta] Sum[(\[Omega]0^2+Abs[2Pi n/(\[Beta]\[Beta] hbar)]\[Gamma]hat[Abs[2Pi n/(\[Beta]\[Beta] hbar)]])/(\[Omega]0^2+(2Pi n/(\[Beta]\[Beta] hbar))^2+Abs[2Pi n/(\[Beta]\[Beta] hbar)]\[Gamma]hat[Abs[2Pi n/(\[Beta]\[Beta] hbar)]]),{n,-1000,1000}];
GLaplace[z_]:=1/(\[Omega]0^2+z^2+z \[Gamma]hat[z]);
GFourier[\[Omega]_]:=GLaplace[-I \[Omega]]-GLaplace[I \[Omega]];
(*Funciones del ba\[NTilde]o*)
SFourier[\[Omega]_]:=I hbar/(2M)Coth[hbar \[Omega] \[Beta]\[Beta]/2]GFourier[\[Omega]];
SDesdeFourierNumerico[t_]:=-NIntegrate[SFourier[\[Omega]]Cos[\[Omega] t],{\[Omega],-Infinito,Infinito},Method->{"LocalAdaptive","SingularityDepth"->20,"SingularityHandler"->Automatic}]/(2 Pi)//Chop;
SDesdeFourierAnalitico[\[Tau]_]:=-I /(M \[Beta]\[Beta] Sqrt[2 Pi])ParallelSum[Chop[term/.n->i/.t->\[Tau]],{i,-sum,sum},DistributedContexts->All]//Chop;
ADesdeLaplace[t_]:=(-hbar/(2 M))InverseLaplaceTransform[GLaplace[z],z,\[Tau]]/.\[Tau]->t//Chop;
ADesdeFourier[t_]:=(-hbar/(2 M Sqrt[2 Pi]))InverseFourierTransform[GFourier[\[Omega]],\[Omega],\[Tau]]/.\[Tau]->t//Chop;
A[t_]:=ADesdeLaplace[t];
S[t_]:=SDesdeFourierNumerico[t];
DA[t_]:=(A[t+0.5*mesh]-A[t-0.5*mesh])/mesh;
DS[t_]:=(S[t+0.5*mesh]-S[t-0.5*mesh])/mesh;
DDA[t_]:=(A[t+mesh]-2.0* A[t]+A[t-mesh])/mesh^2;
DDS[t_]:=(S[t+mesh]-2.0* S[t]+S[t-mesh])/mesh^2;
\[Gamma]Drude[\[Omega]_]:=\[Gamma]\[Gamma] \[Omega]D/(\[Omega]+\[Omega]D);


CompositionLaw[A_,B_]:={A[[1]].B[[1]],A[[1]].B[[2]].Transpose[A[[1]]]+A[[2]]};
InverseTN[A_]:=Module[{a},
a=Inverse[A[[1]]]//Chop;
{a,-a.A[[2]].Transpose[a]}
];
cptpCondition[A_]:=PositiveSemidefiniteMatrixCustom2Q[cptpMatrix[A]];
cptpMatrix[A_]:=Chop[A[[2]]-0.5*I \[CapitalOmega]+0.5*I A[[1]].\[CapitalOmega].Transpose[A[[1]]]];
ComputeF[A_]:=Module[{list},
list=Chop[Eigenvalues[cptpMatrix[A]]];
0.5*Total[Re[Abs[list]-list]]
];
(*EvaluateFunctionCollection[A_,DA_,DDA_,S_,DS_,DDS_,qsq_,psq_]:=Module[{\[Alpha],\[Beta],\[Theta],\[Gamma],\[Delta],\[Eta]},
\[Alpha]=M DA/A;
\[Beta]=hbar/(2 A);
\[Theta]=2/hbar M^2(DDA-DA^2/A);
\[Gamma]=psq/(2 hbar)-M DS/(2 A)+hbar qsq/(8A^2)(1-S^2/qsq^2);
\[Delta]=-M^2/hbar(DS DA/A-DDS)+M*qsq/(2A^2)(DA(1-S^2/qsq^2)+A S DS/qsq^2);
\[Eta]= M^2*DA^2/(2 hbar A^2)*qsq+psq/(2 hbar)-M^2/(2 hbar qsq)(DA/A S-DS)^2;
{\[Alpha],\[Beta],\[Theta],\[Gamma],\[Delta],\[Eta]}
];*)
legends={0,"A","DA","DDA","S","DS","DDS"};
legendscomponents={"T11","T12","T21","T22","N11","N12","N21","N22"};
PlottingCorrelations[listcorr_,which_]:=Module[{ones,target,len},
len=Length[which];
ones=ConstantArray[1,len];
target=Riffle[ones,which];
Transpose[Map[Partition[#,2]&,Transpose[listcorr][[All,target]]]]
];
PlottingComponents[listcorr_,which_]:=Transpose[Map[First[{Partition[Riffle[ConstantArray[#[[1]],Length[which]],Flatten[#[[{2,3}]]][[which]]],2]}]&,ComputeTN[listcorr]]];
PlottingComponentsInverse[listcorr_,which_]:=Transpose[Map[First[{Partition[Riffle[ConstantArray[#[[1]],Length[which]],Flatten[#[[2]]][[which]]],2]}]&,ComputeTNInverse[listcorr]]];

(*Contstruye una rutina q de putazo ya haga los plots para el set de parametros.*)

ComputeFlist[listcorr_]:=Module[{list},
list=ComputeTN[listcorr];
Table[{list[[i]][[1]],ComputeF[CompositionLaw[list[[i+1]][[{2,3}]],InverseTN[list[[i]][[{2,3}]]]//Chop]//Chop]},{i,1,Length[list]-1}]
];
ComputecptpMatrixlist[listcorr_]:=Module[{list},
list=ComputeTN[listcorr];
Table[{list[[i]][[1]],cptpMatrix[CompositionLaw[list[[i+1]][[{2,3}]],InverseTN[list[[i]][[{2,3}]]]//Chop]//Chop]//Chop},{i,1,Length[list]-1}]
];
ComputecptpMatrixComponents[listcorr_,which_]:=Transpose[Map[First[{Partition[Riffle[ConstantArray[#[[1]],Length[which]],Flatten[#[[{2}]]][[which]]],2]}]&,ComputecptpMatrixlist[listcorr]]];
ComputeTN[list_]:=Module[{q,p},
q=qsq;
p=psq;
Map[Chop[{#[[1]],matrixTstandAlone@@#[[{2,3,4}]],matrixNstandAlone@@Flatten[{Take[#,{2,7}],{q,p}}]}]&,Transpose[ list]]];
ComputeTNInverse[list_]:=Map[{#[[1]],InverseTN[#[[{2,3}]]]}&,ComputeTN[list]];



(* ::Subsubsubsection::Closed:: *)
(*term Drude*)


(* ::Input::Initialization:: *)
term:=(I E^(-I t Root[I \[Omega]0^2 \[Omega]D+(-\[Omega]0^2-\[Gamma]\[Gamma] \[Omega]D) #1-I \[Omega]D #1^2+#1^3&,1]) hbar^2 Sqrt[\[Pi]/2] \[Beta]\[Beta]^2 Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,1] (-I \[Omega]D+Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,1]) (Sign[t] (-1+Sign[Abs[Im[Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,1]]]])+2 HeavisideTheta[-t Sign[Im[Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,1]]]] Sign[Im[Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,1]]]))/((-2 I n \[Pi]+hbar \[Beta]\[Beta] Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,1]) (2 I n \[Pi]+hbar \[Beta]\[Beta] Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,1]) (Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,1]-Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,2]) (Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,1]-Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,3]))+(I E^(-I t Root[I \[Omega]0^2 \[Omega]D+(-\[Omega]0^2-\[Gamma]\[Gamma] \[Omega]D) #1-I \[Omega]D #1^2+#1^3&,2]) hbar^2 Sqrt[\[Pi]/2] \[Beta]\[Beta]^2 Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,2] (-I \[Omega]D+Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,2]) (Sign[t] (-1+Sign[Abs[Im[Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,2]]]])+2 HeavisideTheta[-t Sign[Im[Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,2]]]] Sign[Im[Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,2]]]))/((-Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,1]+Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,2]) (-2 I n \[Pi]+hbar \[Beta]\[Beta] Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,2]) (2 I n \[Pi]+hbar \[Beta]\[Beta] Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,2]) (Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,2]-Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,3]))+(I E^(-I t Root[I \[Omega]0^2 \[Omega]D+(-\[Omega]0^2-\[Gamma]\[Gamma] \[Omega]D) #1-I \[Omega]D #1^2+#1^3&,3]) hbar^2 Sqrt[\[Pi]/2] \[Beta]\[Beta]^2 Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,3] (-I \[Omega]D+Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,3]) (Sign[t] (-1+Sign[Abs[Im[Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,3]]]])+2 HeavisideTheta[-t Sign[Im[Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,3]]]] Sign[Im[Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,3]]]))/((-Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,1]+Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,3]) (-Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,2]+Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,3]) (-2 I n \[Pi]+hbar \[Beta]\[Beta] Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,3]) (2 I n \[Pi]+hbar \[Beta]\[Beta] Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,3]))-(I E^(-I t Root[-I \[Omega]0^2 \[Omega]D+(-\[Omega]0^2-\[Gamma]\[Gamma] \[Omega]D) #1+I \[Omega]D #1^2+#1^3&,1]) hbar^2 Sqrt[\[Pi]/2] \[Beta]\[Beta]^2 Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,1] (I \[Omega]D+Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,1]) (Sign[t] (-1+Sign[Abs[Im[Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,1]]]])+2 HeavisideTheta[-t Sign[Im[Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,1]]]] Sign[Im[Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,1]]]))/((-2 I n \[Pi]+hbar \[Beta]\[Beta] Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,1]) (2 I n \[Pi]+hbar \[Beta]\[Beta] Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,1]) (Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,1]-Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,2]) (Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,1]-Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,3]))-(I E^(-I t Root[-I \[Omega]0^2 \[Omega]D+(-\[Omega]0^2-\[Gamma]\[Gamma] \[Omega]D) #1+I \[Omega]D #1^2+#1^3&,2]) hbar^2 Sqrt[\[Pi]/2] \[Beta]\[Beta]^2 Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,2] (I \[Omega]D+Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,2]) (Sign[t] (-1+Sign[Abs[Im[Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,2]]]])+2 HeavisideTheta[-t Sign[Im[Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,2]]]] Sign[Im[Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,2]]]))/((-Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,1]+Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,2]) (-2 I n \[Pi]+hbar \[Beta]\[Beta] Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,2]) (2 I n \[Pi]+hbar \[Beta]\[Beta] Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,2]) (Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,2]-Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,3]))-(I E^(-I t Root[-I \[Omega]0^2 \[Omega]D+(-\[Omega]0^2-\[Gamma]\[Gamma] \[Omega]D) #1+I \[Omega]D #1^2+#1^3&,3]) hbar^2 Sqrt[\[Pi]/2] \[Beta]\[Beta]^2 Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,3] (I \[Omega]D+Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,3]) (Sign[t] (-1+Sign[Abs[Im[Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,3]]]])+2 HeavisideTheta[-t Sign[Im[Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,3]]]] Sign[Im[Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,3]]]))/((-Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,1]+Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,3]) (-Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,2]+Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,3]) (-2 I n \[Pi]+hbar \[Beta]\[Beta] Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,3]) (2 I n \[Pi]+hbar \[Beta]\[Beta] Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,3]))-(E^((2 n \[Pi] t)/(hbar \[Beta]\[Beta])) hbar^2 Sqrt[\[Pi]/2] \[Beta]\[Beta]^2 (2 n \[Pi]-hbar \[Beta]\[Beta] \[Omega]D) (Sign[t] (-1+Sign[Abs[Re[n/(hbar \[Beta]\[Beta])]]])+2 HeavisideTheta[-t Sign[Re[n/(hbar \[Beta]\[Beta])]]] Sign[Re[n/(hbar \[Beta]\[Beta])]]))/(2 (2 n \[Pi]+I hbar \[Beta]\[Beta] Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,1]) (-2 I n \[Pi]+hbar \[Beta]\[Beta] Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,2]) (2 n \[Pi]+I hbar \[Beta]\[Beta] Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,3]))+(E^((2 n \[Pi] t)/(hbar \[Beta]\[Beta])) hbar^2 Sqrt[\[Pi]/2] \[Beta]\[Beta]^2 (2 n \[Pi]+hbar \[Beta]\[Beta] \[Omega]D) (Sign[t] (-1+Sign[Abs[Re[n/(hbar \[Beta]\[Beta])]]])+2 HeavisideTheta[-t Sign[Re[n/(hbar \[Beta]\[Beta])]]] Sign[Re[n/(hbar \[Beta]\[Beta])]]))/(2 (2 n \[Pi]+I hbar \[Beta]\[Beta] Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,1]) (-2 I n \[Pi]+hbar \[Beta]\[Beta] Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,2]) (2 n \[Pi]+I hbar \[Beta]\[Beta] Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,3]))+(E^(-((2 n \[Pi] t)/(hbar \[Beta]\[Beta]))) hbar^2 Sqrt[\[Pi]/2] \[Beta]\[Beta]^2 (2 n \[Pi]+hbar \[Beta]\[Beta] \[Omega]D) (Sign[t] (-1+Sign[Abs[Re[n/(hbar \[Beta]\[Beta])]]])-2 HeavisideTheta[t Sign[Re[n/(hbar \[Beta]\[Beta])]]] Sign[Re[n/(hbar \[Beta]\[Beta])]]))/(2 (2 n \[Pi]-I hbar \[Beta]\[Beta] Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,1]) (2 I n \[Pi]+hbar \[Beta]\[Beta] Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,2]) (2 n \[Pi]-I hbar \[Beta]\[Beta] Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,3]))-(E^(-((2 n \[Pi] t)/(hbar \[Beta]\[Beta]))) hbar^2 Sqrt[\[Pi]/2] \[Beta]\[Beta]^2 (2 n \[Pi]-hbar \[Beta]\[Beta] \[Omega]D) (Sign[t] (-1+Sign[Abs[Re[n/(hbar \[Beta]\[Beta])]]])-2 HeavisideTheta[t Sign[Re[n/(hbar \[Beta]\[Beta])]]] Sign[Re[n/(hbar \[Beta]\[Beta])]]))/(2 (2 n \[Pi]-I hbar \[Beta]\[Beta] Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,1]) (2 I n \[Pi]+hbar \[Beta]\[Beta] Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,2]) (2 n \[Pi]-I hbar \[Beta]\[Beta] Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,3]))


(* ::Subsection:: *)
(*Rest of the package*)


End[]
EndPackage[]
