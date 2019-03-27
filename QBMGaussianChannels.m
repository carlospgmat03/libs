(* ::Package:: *)

BeginPackage["QBMGaussianChannels`",{"Quantum`","Quantum2`"}]
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
PlotA::usage = "PlotA."
PlotS::usage = "PlotS."
PlotN::usage = "PlotN."
PlotT::usage= "PlotT."
PlotNSingular::usage = "Plots singular values of N"
PlotTSingular::usage= "Plots singular values of T"
ComputeSingularValues::usage = "ComputeSingularValues[listcorr_,which_]"
PlotTN::usage = "PlotTN."
PlotNInverse::usage = "PlotNInverse."
PlotTInverse::usage= "PlotTInverse."
PlotTNInverse::usage = "PlotTNInverse."
PlotF::usage = "PlotF."
Correlatorold::usage = "Correaltor[init_,end_,step_] Constructs listcorr."
Correlator::usage = "Correlator[init_,end_,step_,o_:1,d_:4]."
BarridoEn\[Omega]0y\[Gamma]\[Gamma]::usage = "BarridoEn\[Omega]0y\[Gamma]\[Gamma][limite\[Gamma]1_,limite\[Gamma]2_,delta\[Gamma]_,limite\[Omega]01_,limite\[Omega]02_,delta\[Omega]0_,limitetiempo1_,limitetiempo2_,step_]."
DrudeTable::usage = "DrudeTable[init_,end_,step_,time_]."
CustomPseudoInverse::usage = "CustomInverse[m_,tolerance_]"
ComputeTN\[Epsilon]list::usage = "ComputeTN\[Epsilon]list[listcorr_]"
ComputeFHandlingAsymptoticLimit::usage = "listFHandlingAsymptoticLimit[listcorr_,epsilon_:0.0001]"
\[CapitalLambda]::usage = "Matrix to get rid of the units."
TimeScaleCorrector::usage = "TimeScaleCorrector[factor_,list_]."


(* ::Input::Initialization:: *)
TimeScaleCorrector[factor_,paso_,list_]:=Module[{},({#1[[1]]paso,#1[[2]]factor}&)/@list];
\[CapitalOmega]={{0,1},{-1,0}};
matrixTstandAlone[A_,DA_,DDA_]:={{-((2 M DA)/hbar),(2 A)/hbar},{(2 M^2 DDA)/hbar,-((2 M DA)/hbar)}}//Chop;
matrixNstandAlone[A_,DA_,DDA_,S_,DS_,DDS_,qsq_,psq_]:={{qsq+(4 A (psq A-hbar M DS))/hbar^2-S^2/qsq,-((2 M (A (-2 psq DA+hbar M DDS)+hbar M DA DS))/hbar^2)-(M DS S)/qsq},{-((2 M (A (-2 psq DA+hbar M DDS)+hbar M DA DS))/hbar^2)-(M DS S)/qsq,psq+M^2 ((4 DA (psq DA-hbar M DDS))/hbar^2-DS^2/qsq)}}//Chop;
matrixTstandAlone[t_]:={{-((2 M DA[t])/hbar),(2 A[t])/hbar},{(2 M^2 DDA[t])/hbar,-((2 M DA[t])/hbar)}};
matrixNstandAlone[t_]:={{qsq+(4 A[t] (psq A[t]-hbar M DS[t]))/hbar^2-S[t]^2/qsq,-((2 M (A[t] (-2 psq DA[t]+hbar M DDS[t])+hbar M DA[t] DS[t]))/hbar^2)-(M DS[t] S[t])/qsq},{-((2 M (A[t] (-2 psq DA[t]+hbar M DDS[t])+hbar M DA[t] DS[t]))/hbar^2)-(M DS[t] S[t])/qsq,psq+M^2 ((4 DA[t] (psq DA[t]-hbar M DDS[t]))/hbar^2-DS[t]^2/qsq)}};
\[HBar]:=hbar;

(*Desviaciones*)
qsq:=1/(M \[Beta]\[Beta]) Sum[1/(\[Omega]0^2+(2Pi n/(\[Beta]\[Beta] hbar))^2+Abs[2Pi n/(\[Beta]\[Beta] hbar)] \[Gamma]hat[Abs[2Pi n/(\[Beta]\[Beta] hbar)]]),{n,-10000,1000}];
psq:=M/\[Beta]\[Beta] Sum[(\[Omega]0^2+Abs[2Pi n/(\[Beta]\[Beta] hbar)]\[Gamma]hat[Abs[2Pi n/(\[Beta]\[Beta] hbar)]])/(\[Omega]0^2+(2Pi n/(\[Beta]\[Beta] hbar))^2+Abs[2Pi n/(\[Beta]\[Beta] hbar)]\[Gamma]hat[Abs[2Pi n/(\[Beta]\[Beta] hbar)]]),{n,-10000,10000}];
GLaplace[z_]:=1/(\[Omega]0^2+z^2+z \[Gamma]hat[z]);
GFourier[\[Omega]_]:=GLaplace[-I \[Omega]]-GLaplace[I \[Omega]];
(*Funciones del ba\[NTilde]o*)
SFourier[\[Omega]_]:=I hbar/(2M)Coth[hbar \[Omega] \[Beta]\[Beta]/2]GFourier[\[Omega]];
SDesdeFourierNumerico[t_]:=-NIntegrate[SFourier[\[Omega]]Cos[\[Omega] t],{\[Omega],-Infinito,Infinito},Method->{"LevinRule","Kernel"->Cos[\[Omega] t]}(*,MaxRecursion->20*)]/(2 Pi)//Chop;
SDesdeFourierAnalitico[\[Tau]_]:=-I /(M \[Beta]\[Beta] Sqrt[2 Pi])Sum[Chop[Chop[term]/.n->i/.t->\[Tau]],{i,-sum,sum}]//Chop;
ADesdeLaplace[t_]:=(-hbar/(2 M))InverseLaplaceTransform[GLaplace[z],z,\[Tau]]/.\[Tau]->t//Chop;
ADesdeFourier[t_]:=(-hbar/(2 M Sqrt[2 Pi]))InverseFourierTransform[GFourier[\[Omega]],\[Omega],\[Tau]]/.\[Tau]->t//Chop;
A[t_]:=ADesdeLaplace[t];
S[t_]:=SDesdeFourierNumerico[t];
DA[t_]:=(A[t+0.5*mesh]-A[t-0.5*mesh])/mesh;
DS[t_]:=(S[t+0.5*mesh]-S[t-0.5*mesh])/mesh;
DDA[t_]:=(A[t+mesh]-2.0* A[t]+A[t-mesh])/mesh^2;
DDS[t_]:=(S[t+mesh]-2.0* S[t]+S[t-mesh])/mesh^2;
\[Gamma]Drude[\[Omega]_]:=\[Gamma]\[Gamma] \[Omega]D/(\[Omega]+\[Omega]D);


CompositionLaw[A_,B_]:=Chop[{A[[1]].B[[1]],A[[1]].B[[2]].Transpose[A[[1]]]+A[[2]]}];
InverseTN[A_]:=Module[{a,d},
d=0.00000001;
a=CustomPseudoInverse[A[[1]],d];
{a[[1]],If[a[[2]],{{0,0},{0,0}},-a[[1]].A[[2]].Transpose[a[[1]]]]}
(*a=PseudoInverse[A[[1]],0.001];
{a,-a.A[[2]].Transpose[a]}*)
]
cptpCondition[A_]:=PositiveSemidefiniteMatrixCustom2Q[cptpMatrix[A]];
cptpMatrix[A_]:=Chop[\[CapitalLambda].A[[2]].\[CapitalLambda]-0.5*I*\[CapitalOmega]+ hbar 0.5*I*\[CapitalLambda].A[[1]].\[CapitalOmega].Transpose[A[[1]]].\[CapitalLambda]];
ComputeF[A_]:=Module[{list},
list=Re[Chop[Eigenvalues[cptpMatrix[A]]]];
0.5*Total[Abs[list]-list]
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

PlottingSingular[listcorr_,which_]:=Module[{list,tmp},
list=ComputeTN[listcorr];
Transpose[Table[tmp=SingularValueList[list[[i]][[which+1]]];{{list[[i]][[1]],tmp[[1]]},{list[[i]][[1]],tmp[[2]]}},{i,1,Length[list]}]]
];(*Map[Partition[Riffle[ConstantArray[#[[1]],2],SingularValueList[#[[1+which]]]],2]&,ComputeTN[listcorr]];*)

ComputeSingularValues[listcorr_,which_]:=Module[{list,tmp},
list=ComputeTN[listcorr];
Transpose[Table[tmp=SingularValueList[list[[i]][[which+1]]];{{list[[i]][[1]],tmp[[1]]},{list[[i]][[1]],tmp[[2]]}},{i,1,Length[list]}]]
];

PlotTSingular:=ListLinePlot[PlottingSingular[listcorr,1]];

PlotNSingular:=ListLinePlot[PlottingSingular[listcorr,2]];

ComputeFHandlingAsymptoticLimit[listcorr_,epsilon_:0.0001]:=Module[{out,list,list1,list2,pos1,pos2,limit,listF},
list=ComputeSingularValues[listcorr,1];
listF=ComputeFlist[listcorr];
list1=list[[1]];
list2=list[[2]];
pos1=Position[list1[[All,1]],Select[list1,Abs[#[[2]]]<epsilon&][[1,1]]][[1,1]];
pos2=Position[list2[[All,1]],Select[list2,Abs[#[[2]]]<epsilon&][[1,1]]][[1,1]];
limit=Min[{pos1,pos2}];
Table[listF[[i]]=Flatten[{listF[[i]][[1]],listF[[limit]][[2]]}],{i,limit,Length[listF]}];
listF
]

ComputeFlist[listcorr_]:=Module[{list},
list=Chop[ComputeTN[listcorr]];
Table[{list[[i]][[1]],ComputeF[CompositionLaw[list[[i+1]][[{2,3}]],InverseTN[list[[i]][[{2,3}]]]//Chop]//Chop]},{i,1,Length[list]-1}]
];
ComputeTN\[Epsilon]list[listcorr_]:=Module[{list},
list=Chop[ComputeTN[listcorr]];
Table[{list[[i]][[1]],CompositionLaw[list[[i+1]][[{2,3}]],InverseTN[list[[i]][[{2,3}]]]//Chop]},{i,1,Length[list]-1}]
];
ComputecptpMatrixlist[listcorr_]:=Module[{list},
list=ComputeTN[listcorr];
Table[{list[[i]][[1]],cptpMatrix[CompositionLaw[list[[i+1]][[{2,3}]],InverseTN[list[[i]][[{2,3}]]]//Chop]//Chop]//Chop},{i,1,Length[list]-1}]
];
ComputecptpMatrixComponents[listcorr_,which_]:=Transpose[Map[First[{Partition[Riffle[ConstantArray[#[[1]],Length[which]],Flatten[#[[{2}]]][[which]]],2]}]&,ComputecptpMatrixlist[listcorr]]];
ComputeTN[list_]:=Module[{q,p},
q=qsq;
p=psq;
Map[Chop[{#[[1]],Threshold[matrixTstandAlone@@#[[{2,3,4}]],0.0000001],Chop[matrixNstandAlone@@Flatten[{Take[#,{2,7}],{q,p}}],0.00001]}]&,Transpose[ list]]];

ComputeTNInverse[list_]:=Map[{#[[1]],InverseTN[#[[{2,3}]]]}&,ComputeTN[list]];
PlotA:=ListLinePlot[PlottingCorrelations[listcorr,pl={2,3,4}],PlotRange->All,PlotLegends->Placed[legends[[pl]],Above]];
PlotS:=ListLinePlot[PlottingCorrelations[listcorr,pl={5,6,7}],PlotRange->All,PlotLegends->Placed[legends[[pl]],Above]];
PlotN:=ListLinePlot[PlottingComponents[listcorr,pl={5,6,7,8}],PlotRange->All,PlotLegends->Placed[legendscomponents[[pl]],Right],PlotStyle->Automatic];
PlotTN:=ListLinePlot[PlottingComponents[listcorr,pl={1,2,3,4,5,6,7,8}],PlotRange->Automatic,PlotLegends->Placed[legendscomponents[[pl]],Right],PlotStyle->{Automatic,Automatic,Automatic,Automatic,Automatic,Automatic,Automatic,Directive[Dashed,Opacity[0.5]]}];
PlotT:=ListLinePlot[PlottingComponents[listcorr,pl={1,2,3,4}],PlotRange->All,PlotLegends->Placed[legendscomponents[[pl]],Right],PlotStyle->{Automatic,Automatic,Automatic,Automatic,Automatic,Automatic,Automatic,Directive[Dashed,Opacity[0.5]]}];
PlotNInverse:=ListLinePlot[PlottingComponentsInverse[listcorr,pl={5,6,7,8}],PlotRange->All,PlotLegends->Placed[legendscomponents[[pl]],Right],PlotStyle->Automatic];
PlotTNInverse:=ListLinePlot[PlottingComponentsInverse[listcorr,pl={1,2,3,4,5,6,7,8}],PlotRange->Automatic,PlotLegends->Placed[legendscomponents[[pl]],Right],PlotStyle->{Automatic,Automatic,Automatic,Automatic,Automatic,Automatic,Automatic,Directive[Dashed,Opacity[0.5]]}];
PlotTInverse:=ListLinePlot[PlottingComponentsInverse[listcorr,pl={1,2,3,4}],PlotRange->All,PlotLegends->Placed[legendscomponents[[pl]],Right],PlotStyle->{Automatic,Automatic,Automatic,Automatic,Automatic,Automatic,Automatic,Directive[Dashed,Opacity[0.5]]}];
PlotF:=ListLinePlot[listF,PlotRange->All,PlotStyle->Red];

Correlator[init_,end_,step_,o_:1,d_:4]:=Module[{list,h,x,DA,DDA,DS,DDS},
list=Table[{i,A[i],S[i]},{i,init,end+2*o*step,step}];
DA=GaussianFilter[ForwardNumericalDifferenceForTimeSeries[list[[All,2]],o]/step,d];
DDA=GaussianFilter[ForwardNumericalDifferenceForTimeSeries[DA,o]/step,d];
DS=GaussianFilter[ForwardNumericalDifferenceForTimeSeries[list[[All,3]],o]/step,d];
DDS=GaussianFilter[ForwardNumericalDifferenceForTimeSeries[DS,o]/step,d];
DA=Take[DA,Length[list]-2*o];
DS=Take[DS,Length[list]-2*o];
list=Take[list,Length[list]-2*o];
list=Transpose[list];
{list[[1]],list[[2]],DA,DDA,list[[3]],DS,DDS}
];

Correlatorold[init_,end_,step_]:=Module[{},
listcorr=Table[{i,A[i],DA[i],DDA[i],S[i],DS[i],DDS[i]},{i,init,end,step}]//Transpose];

BarridoEn\[Omega]0y\[Gamma]\[Gamma][limite\[Gamma]1_,limite\[Gamma]2_,delta\[Gamma]_,limite\[Omega]01_,limite\[Omega]02_,delta\[Omega]0_,limitetiempo1_,limitetiempo2_,step_]:=Module[{F,listcorr,FforCharlie,Flist},
Flatten[ParallelTable[
listcorr=Table[{i,A[i],DA[i],DDA[i],S[i],DS[i],DDS[i]},{i,limitetiempo1,limitetiempo2,step}]//Transpose;
Flist=ComputeFlist[listcorr];
F=Flist[[All,2]];
FforCharlie=ListIntegrate[Flist][[All,2]];
{\[Gamma]\[Gamma]=i,\[Omega]0=j,Transpose[listcorr],step*Total[F],CharlieMeasure[FforCharlie]}
,{i,limite\[Gamma]1,limite\[Gamma]2,delta\[Gamma]},{j,limite\[Omega]01,limite\[Omega]02,delta\[Omega]0},DistributedContexts->All]
,1]];
DrudeTable[init_,end_,step_,timeinit_,timeend_,timestep_]:=ParallelTable[\[Omega]D=j;ComputeFHandlingAsymptoticLimit[Table[{i,A[i],DA[i],DDA[i],S[i],DS[i],DDS[i]},{i,timeinit,timeend,timestep}]//Transpose],{j,init,end,step},DistributedContexts->All];
DrudeTable[init_,end_,step_]:=DrudeTable[init,end,step,0.01,20.0,0.1];
CustomPseudoInverse[m_,tolerance_]:=Module[{s,v,d,vinv,f},
{s,v,d}=SingularValueDecomposition[m];
f=False;
vinv=SparseArray[{i_, i_}:>If[v[[i,i]]>tolerance,1/v[[i,i]],f=True;0],Length[v]];
{d.vinv.ConjugateTranspose[s],f}
];
(*Para quitar unidades*)
q0=Sqrt[hbar/(M \[Omega]0)];p0=Sqrt[hbar M \[Omega]0];
\[CapitalLambda]:=DiagonalMatrix[{q0,p0}]//Inverse;


(* ::Subsubsubsection::Closed:: *)
(*term Drude*)


(* ::Input::Initialization:: *)
term:=(I E^(-I t Root[I \[Omega]0^2 \[Omega]D+(-\[Omega]0^2-\[Gamma]\[Gamma] \[Omega]D) #1-I \[Omega]D #1^2+#1^3&,1]) hbar^2 Sqrt[\[Pi]/2] \[Beta]\[Beta]^2 Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,1] (-I \[Omega]D+Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,1]) (Sign[t] (-1+Sign[Abs[Im[Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,1]]]])+2 HeavisideTheta[-t Sign[Im[Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,1]]]] Sign[Im[Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,1]]]))/((-2 I n \[Pi]+hbar \[Beta]\[Beta] Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,1]) (2 I n \[Pi]+hbar \[Beta]\[Beta] Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,1]) (Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,1]-Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,2]) (Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,1]-Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,3]))+(I E^(-I t Root[I \[Omega]0^2 \[Omega]D+(-\[Omega]0^2-\[Gamma]\[Gamma] \[Omega]D) #1-I \[Omega]D #1^2+#1^3&,2]) hbar^2 Sqrt[\[Pi]/2] \[Beta]\[Beta]^2 Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,2] (-I \[Omega]D+Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,2]) (Sign[t] (-1+Sign[Abs[Im[Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,2]]]])+2 HeavisideTheta[-t Sign[Im[Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,2]]]] Sign[Im[Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,2]]]))/((-Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,1]+Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,2]) (-2 I n \[Pi]+hbar \[Beta]\[Beta] Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,2]) (2 I n \[Pi]+hbar \[Beta]\[Beta] Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,2]) (Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,2]-Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,3]))+(I E^(-I t Root[I \[Omega]0^2 \[Omega]D+(-\[Omega]0^2-\[Gamma]\[Gamma] \[Omega]D) #1-I \[Omega]D #1^2+#1^3&,3]) hbar^2 Sqrt[\[Pi]/2] \[Beta]\[Beta]^2 Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,3] (-I \[Omega]D+Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,3]) (Sign[t] (-1+Sign[Abs[Im[Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,3]]]])+2 HeavisideTheta[-t Sign[Im[Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,3]]]] Sign[Im[Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,3]]]))/((-Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,1]+Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,3]) (-Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,2]+Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,3]) (-2 I n \[Pi]+hbar \[Beta]\[Beta] Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,3]) (2 I n \[Pi]+hbar \[Beta]\[Beta] Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,3]))-(I E^(-I t Root[-I \[Omega]0^2 \[Omega]D+(-\[Omega]0^2-\[Gamma]\[Gamma] \[Omega]D) #1+I \[Omega]D #1^2+#1^3&,1]) hbar^2 Sqrt[\[Pi]/2] \[Beta]\[Beta]^2 Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,1] (I \[Omega]D+Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,1]) (Sign[t] (-1+Sign[Abs[Im[Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,1]]]])+2 HeavisideTheta[-t Sign[Im[Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,1]]]] Sign[Im[Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,1]]]))/((-2 I n \[Pi]+hbar \[Beta]\[Beta] Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,1]) (2 I n \[Pi]+hbar \[Beta]\[Beta] Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,1]) (Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,1]-Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,2]) (Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,1]-Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,3]))-(I E^(-I t Root[-I \[Omega]0^2 \[Omega]D+(-\[Omega]0^2-\[Gamma]\[Gamma] \[Omega]D) #1+I \[Omega]D #1^2+#1^3&,2]) hbar^2 Sqrt[\[Pi]/2] \[Beta]\[Beta]^2 Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,2] (I \[Omega]D+Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,2]) (Sign[t] (-1+Sign[Abs[Im[Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,2]]]])+2 HeavisideTheta[-t Sign[Im[Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,2]]]] Sign[Im[Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,2]]]))/((-Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,1]+Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,2]) (-2 I n \[Pi]+hbar \[Beta]\[Beta] Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,2]) (2 I n \[Pi]+hbar \[Beta]\[Beta] Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,2]) (Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,2]-Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,3]))-(I E^(-I t Root[-I \[Omega]0^2 \[Omega]D+(-\[Omega]0^2-\[Gamma]\[Gamma] \[Omega]D) #1+I \[Omega]D #1^2+#1^3&,3]) hbar^2 Sqrt[\[Pi]/2] \[Beta]\[Beta]^2 Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,3] (I \[Omega]D+Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,3]) (Sign[t] (-1+Sign[Abs[Im[Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,3]]]])+2 HeavisideTheta[-t Sign[Im[Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,3]]]] Sign[Im[Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,3]]]))/((-Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,1]+Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,3]) (-Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,2]+Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,3]) (-2 I n \[Pi]+hbar \[Beta]\[Beta] Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,3]) (2 I n \[Pi]+hbar \[Beta]\[Beta] Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,3]))-(E^((2 n \[Pi] t)/(hbar \[Beta]\[Beta])) hbar^2 Sqrt[\[Pi]/2] \[Beta]\[Beta]^2 (2 n \[Pi]-hbar \[Beta]\[Beta] \[Omega]D) (Sign[t] (-1+Sign[Abs[Re[n/(hbar \[Beta]\[Beta])]]])+2 HeavisideTheta[-t Sign[Re[n/(hbar \[Beta]\[Beta])]]] Sign[Re[n/(hbar \[Beta]\[Beta])]]))/(2 (2 n \[Pi]+I hbar \[Beta]\[Beta] Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,1]) (-2 I n \[Pi]+hbar \[Beta]\[Beta] Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,2]) (2 n \[Pi]+I hbar \[Beta]\[Beta] Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,3]))+(E^((2 n \[Pi] t)/(hbar \[Beta]\[Beta])) hbar^2 Sqrt[\[Pi]/2] \[Beta]\[Beta]^2 (2 n \[Pi]+hbar \[Beta]\[Beta] \[Omega]D) (Sign[t] (-1+Sign[Abs[Re[n/(hbar \[Beta]\[Beta])]]])+2 HeavisideTheta[-t Sign[Re[n/(hbar \[Beta]\[Beta])]]] Sign[Re[n/(hbar \[Beta]\[Beta])]]))/(2 (2 n \[Pi]+I hbar \[Beta]\[Beta] Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,1]) (-2 I n \[Pi]+hbar \[Beta]\[Beta] Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,2]) (2 n \[Pi]+I hbar \[Beta]\[Beta] Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,3]))+(E^(-((2 n \[Pi] t)/(hbar \[Beta]\[Beta]))) hbar^2 Sqrt[\[Pi]/2] \[Beta]\[Beta]^2 (2 n \[Pi]+hbar \[Beta]\[Beta] \[Omega]D) (Sign[t] (-1+Sign[Abs[Re[n/(hbar \[Beta]\[Beta])]]])-2 HeavisideTheta[t Sign[Re[n/(hbar \[Beta]\[Beta])]]] Sign[Re[n/(hbar \[Beta]\[Beta])]]))/(2 (2 n \[Pi]-I hbar \[Beta]\[Beta] Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,1]) (2 I n \[Pi]+hbar \[Beta]\[Beta] Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,2]) (2 n \[Pi]-I hbar \[Beta]\[Beta] Root[I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1-I \[Omega]D #1^2+#1^3&,3]))-(E^(-((2 n \[Pi] t)/(hbar \[Beta]\[Beta]))) hbar^2 Sqrt[\[Pi]/2] \[Beta]\[Beta]^2 (2 n \[Pi]-hbar \[Beta]\[Beta] \[Omega]D) (Sign[t] (-1+Sign[Abs[Re[n/(hbar \[Beta]\[Beta])]]])-2 HeavisideTheta[t Sign[Re[n/(hbar \[Beta]\[Beta])]]] Sign[Re[n/(hbar \[Beta]\[Beta])]]))/(2 (2 n \[Pi]-I hbar \[Beta]\[Beta] Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,1]) (2 I n \[Pi]+hbar \[Beta]\[Beta] Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,2]) (2 n \[Pi]-I hbar \[Beta]\[Beta] Root[-I \[Omega]0^2 \[Omega]D-\[Omega]0^2 #1-\[Gamma]\[Gamma] \[Omega]D #1+I \[Omega]D #1^2+#1^3&,3]))


(* ::Subsection:: *)
(*Rest of the package*)


EndPackage[]
