BeginPackage["RankEvolution`"]

Needs["Carlos`"]
(*
Needs["Quantum`"]
*)
(* Model simulation {{{*)
AlgorithStepNormal::usage = "Evolves a system one time step acording to the gaussian scale invariant model. Usage, AlgorithStepNormal[WordOrder_, \[Alpha]_]"
Model::usage = "Gives a full model with arbitrary time steps. Usage, Model[NumberOfWords_, TimeSteps_, alpha_]"
Diversity::usage = "Gives the diversity for a full model. The model is a matrix containing the data per year. Usage: Diversity[data_]"
ProbabilityOfChange::usage = "Gives the probability of change of a Model. Usage: ProbabilityOfChange[data_]"
ListToOrdered::usage = "Gives an unsorted and with any kind of data, a list of nombers that reflects the order"
GenerateMixedModel::usage = "Generates the model that interpolates between the closed and the open one. Usage: GenerateMixedModel[probability_, LengthModel_, EvolutionTime_] or GenerateMixedModel[probability_, LengthModel_, EvolutionTime_, skip_]"	
ChangePosition::usage = "Changes the position according to the null model that is written. Usage, ChangePosition[Psi0_, {a_Integer, b_Integer}], ChangePosition[Psi0_] or ChangePosition[Psi0_, skip_Integer]"
FastProbabilityOfDisplacement::usage = "Calculates the probability of displacement, in a fast way. Usage FastProbabilityOfDisplacement[k_, n_, s_, i_]"
GenerateModel2::usage = "Generates the null model"
ModelStepIonizing::usage = "One step of the ionizing model"
GenerateModelStepIonizing::usage = "Complete set of steps for the ionizing model"
(*}}}*)
Begin["`Private`"]
(* Ionizing Model EL BUENO {{{*)
ModelStepIonizing[ ModelState_List, {ProbabilityLevy_, ProbabilityReplace_}] :=
 Module[{Psi, LatestAddedState},
  {Psi, LatestAddedState} = ModelState;
  Psi = If[Random[] < ProbabilityLevy, ChangePosition[Psi], Psi]; 
  If[Random[] < ProbabilityReplace,
   LatestAddedState = LatestAddedState + 1;
   Psi = ReplacePart[Psi, 
     RandomInteger[{1, Length[Psi]}] -> LatestAddedState];];
  {Psi, LatestAddedState}]

ModelStepIonizing[LengthModel_Integer, ProbabilitySet_List, EvolutionTime_] := 
 Nest[ModelStepIonizing[#, ProbabilitySet] &, {Range[LengthModel], 
    LengthModel}, EvolutionTime][[1]]

ModelStepIonizing[LengthModel_Integer, ProbabilitySet_List] :=
  ModelStepIonizing[LengthModel, ProbabilitySet,1]

GeneratePreModelIonizing[LengthModel_, EvolutionTime_, ProbabilitySet_] :=
 NestList[ModelStepIonizing[#, ProbabilitySet] &, {Range[LengthModel],
    LengthModel}, EvolutionTime - 1]

GenerateModelStepIonizing[LengthModel_, EvolutionTime_, ProbabilitySet_] := 
  Module[{completo},
   completo = GeneratePreModelIonizing[LengthModel, EvolutionTime, ProbabilitySet];
   Transpose[completo][[1]][[All, ;; Length[completo[[1, 1]]]]]];

(* Ionizing Model }}}*)
(* Null Model {{{*)
ChangePosition[Psi0_, {a_Integer, b_Integer}] := Flatten[Which[
   a < b, {Psi0[[1 ;; a - 1]], Psi0[[a + 1 ;; b]], Psi0[[a]], Psi0[[b + 1 ;;]]}, 
   a == b, Psi0, 
   a > b, {Psi0[[1 ;; b - 1]], Psi0[[a]], Psi0[[b ;; a - 1]], Psi0[[a + 1 ;;]]}]]
ChangePosition[Psi0_] := ChangePosition[Psi0, RandomInteger[{1, Length[Psi0]}, {2}]]
ChangePosition[Psi0_, skip_Integer] := Nest[ChangePosition, Psi0, skip]

GenerateModel2[LengthModel_, EvolutionTime_, skip_] := 
 NestList[ChangePosition[#, skip] &, Range[LengthModel], EvolutionTime - 1]
GenerateModel2[LengthModel_, EvolutionTime_] := GenerateModel2[LengthModel, EvolutionTime, 1]

ChangePositionToLast[Psi0_, a_Integer] := ChangePosition[Psi0, {a, Length[Psi0]}]
ChangePositionToLast[Psi0_] := ChangePositionToLast[Psi0, RandomInteger[{1, Length[Psi0]}]]

MixedModelIteration[Psi0_, probability_] := If[Random[] <= probability, ChangePosition[Psi0], ChangePositionToLast[Psi0]]
MixedModelIteration[Psi0_, probability_, skip_Integer] := Nest[MixedModelIteration[#, probability] &, Psi0, skip]
GenerateMixedModel[probability_, LengthModel_, EvolutionTime_, skip_] := NestList[MixedModelIteration[#, probability, skip] &, Range[LengthModel], EvolutionTime - 1]
GenerateMixedModel[probability_, LengthModel_, EvolutionTime_] := GenerateMixedModel[probability, LengthModel, EvolutionTime, 1]

(*}}}*)
(* Model simulation {{{*)
AlgorithStepNormal[WordOrder_, \[Alpha]_] := Module[{NumberOfWords, NewPositions},
  NumberOfWords = Length[WordOrder];
  NewPositions = # + RandomReal[NormalDistribution[0, # \[Alpha]]] & /@ Range[NumberOfWords];
  WordOrder[[Ordering[NewPositions]]]]
Model[NumberOfWords_, TimeSteps_, alpha_] := NestList[AlgorithStepNormal[#, alpha] &, Range[NumberOfWords], TimeSteps]
Unfold[lp2_] := Permute[Range[Length[lp2]], InversePermutation[FindPermutation[lp2, Sort[lp2]]]]
(*}}}*)
(* exact formulae {{{*)
(* En al siguiente:

i: Numero de iteraciones
*)
FastProbabilityOfDisplacement[k_, n_, s_, i_] := If[Abs[s] > i, ProbabilidadLevy[i, n], ConvolutionProbabilityOfDisplacement[k, n, s, i]]

ConvolutionProbabilityOfDisplacement[k_, n_, s_, 1] := ProbabilityOfDisplacement[k, n, s, 1]

ConvolutionProbabilityOfDisplacement[k_, n_, s_, i_] := ConvolutionProbabilityOfDisplacement[k, n, s, i] = 
  N[1/n^2 + Sum[FastProbabilityOfDisplacement[k, n, s - a, i - 1] (ProbabilityOfDisplacement[k + s - a, n, a, 1] - 1./n^2), {a, -1, 1}]]

ProbabilidadLevy[I_, n_] := 1./n (1 - (1 - 1/n)^I)

(*FastProbabilityOfDisplacementCorrectRange[k_, n_, s_, i_] := If[1 <= k <= n && 1 <= (s - k) <= n, FastProbabilityOfDisplacement[k, n, s, i], 0]*)

ProbabilityOfDisplacement[CurrentRankK_, LengthModel_, Displacement_, 1] := ProbabilityOfDisplacement[CurrentRankK, LengthModel, Displacement]

ProbabilityOfDisplacement[CurrentRankK_, LengthModel_, Displacement_] := 
 Levy[CurrentRankK, LengthModel, Displacement] + Drift[CurrentRankK, LengthModel, Displacement]

Levy[k_, n_, s_] := If[1 <= k + s <= n, 1/n^2, 0]
Drift[k_, n_, 1] := (n - k) k/n^2
Drift[k_, n_, -1] := (n - k + 1) (k - 1)/n^2
Drift[k_, n_, 0] := (k - 1)^2/n^2 + (n - k)^2/n^2
Drift[k_, n_, s_] := 0
(*}}}*)
Diversity[data_] := NumberList[N[Length[Union[#]]/Length[#]] & /@ Transpose[data]]
NumberOfChanges[x_List] := Length[Split[x]] - 1
ProbabilityOfChangeSingleData[x_List] := NumberOfChanges[x]/(Length[x] - 1)
ProbabilityOfChange[data_] :=NumberList[ProbabilityOfChangeSingleData /@ Transpose[data]]
ListToOrdered[x_] := x /. (Rule @@ # & /@ Transpose[{Sort[x], Range[Length[x]]}])
End[ ]
EndPackage[ ]

