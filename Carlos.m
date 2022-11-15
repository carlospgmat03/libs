BeginPackage["Carlos`"] ;


ItppVectorToExpression::usage="Helps read data directly from itpp output, like to interface with a cpp program"

MyAverage::usage = 
    "This function gives the average of more complex quantities than lists, 
    for example it is able to process lists of lists";
RandomUnitVector::usage = "This gives a random vector with the Haar measure. The dimension is 
                           the argument. If no argument is suplied, it is asumed to be 3";

DistanceBetweenSetsOfPoints::usage = "It calculates the distance between two sets of points. They might have a different order. Basically it has been used to compare two spectra";

Seed::usage = 
    "This function, without any argument, gives a random integer between 0
    and 1000000000-1 which can be used as a seed for an external program";

Instruction::usage = 
    "This instrucion creats a string that in s Linux or Unix shell will enter
    some given parameters into an external program. The first argument is a list
    of lists with the values to be entered and the second one is the program";

Norma::usage = "This function gives the norm of a list of numbers";

ColorCoding::usage="This recieves 2 integer inputs, and outputs a 
	graphic to read the number that represents each Hue";
(*
Se comenta porque ya Mathematica 13.1 trae una funcion que hace eso
BlockDiagonalMatrix::usage="Gives a block diagonal matrix from a list"
*)

OffSpellErrors::usage ="Turns of spelling errors"
OnSpellErrors::usage ="Turns on spelling errors"

Log10::usage ="Calculates Log10[x_]:=Log[10,x]"
ColumnAddKeepFirst::usage="To add to matrices, keeping the first column of the first matrix untouched"
ReadListUncomment::usage="igual a ReadList[] pero quita todo lo que comienze con #"
NumberList::usage="Number a list, i.e. prepend with an intenger from 1 to the Length of the list"


HistogramListPoints::usage="Shows the points that would correspond to a Histogram. Accepts
the same options as Histogram and HistogramList. Usage HistogramListPoints[data] or
HistogramListPoints[data, bspec] or HistogramList[data,bspec,hspec]"

HistogramPointsForLine::usage="Calculates the points to make a line corresponding to a Histogram. 
Accepts
the same options as Histogram and HistogramList. Usage HistogramPointsForLine [data] or
HistogramPointsForLine[data, bspec] or HistogramPointsForLine[data,bspec,hspec]"



(* {{{ Symbols and legends *)
MySymbol::usage="Para poner simbolos. Tiene defauls. Es el recomendado ahora"
SymbolNumber::usage="Option for MySymbol"
Coordinate::usage="Option for MySymbol"
Color::usage="Option for MySymbol"
Proportion::usage="Option for MySymbol"
delta::usage="Option for MySymbol"
ThicknessBorder::usage="Option for MySymbol"
MyTriangle::usage = "Graphics almost primitive MyTriangle[{x_, y_}, Color1_, Proportion_, delta_, th_] ";
MyInvertedTriangle::usage = "Graphics almost primitive MyInvertedTriangle[{x_, y_}, Color1_, Proportion_, delta_, th_] ";
MySquare::usage = "Graphics almost primitive MySquare[{x_, y_}, Color1_, Proportion_, delta_, th_]";
MyCircle::usage = "Graphics almost primitive MyCircle[{x_, y_}, Color1_, Proportion_, delta_, th_]";
MyRhombous::usage = "Graphics almost primitive MyRhombous[{x_, y_}, Color1_, Proportion_, delta_, th_]";
My4PointStar1::usage = "Graphics almost primitive My4PointStar1[{x_, y_}, Color1_, Proportion_, delta_, th_]";
My4PointStar2::usage = "Graphics almost primitive My4PointStar2[{x_, y_}, Color1_, Proportion_, delta_, th_]";
My4PointStar3::usage = "Graphics almost primitive My4PointStar3[{x_, y_}, Color1_, Proportion_, delta_, th_]";
My5PointStar::usage = "Graphics almost primitive My5PointStar[{x_, y_}, Color1_, Proportion_, delta_, th_]";

InsetWithSymbols::usage="To create nice symbols in plots. "
MyLegend::usage="Ver LegendBox"

LegendBox::usage=" Ejemplos de uso:
{UpperHeight = 0.9, Xpos = .57, Xlength = .1, 
  XSepText = .1, \[CapitalDelta]Height = .15};
kk = LegendBox[{\"b1\", \"b2\", 
    \"b4\"}, {GrayLevel[0], Style\[Beta][#]} & /@ \[Beta]s, UpperHeight,
    Xpos, Xlength, XSepText, \[CapitalDelta]Height];

LegendBox[
 \"n=\" <> ToString[#] & /@ ns, {Thickness[tjh], Hue[#/Length[ns]]} & /@
   Range[Length[ns]], .9, .6, .2, .05, .1]
"

Alignment::usage="Option for LegendBox and MyLegend"
(* }}} *)
(* {{{ Geometry *)
EllipseCharacteristics::usage="Get center, angle of rotation and semiaxis of an elipse. EllipseCharacteristics[poly_, vars_]
 For example, EllipseCharacteristics[4 x^2 - 4 x y + 7 y^2 + 12 x + 6 y - 9, {x, y}]"



(* }}} *)
Begin["Private`"];


(* Read data from itpp output *)
ItppVectorToExpression[vector_String] := 
 ToExpression /@ StringSplit[ StringReplace[
    StringTake[vector, {2, -2}], {"i" -> "I", "e+" :> "*^", "e-" :> "*^-"}]]
(*  Geometry *)
EllipseCharacteristics[poly_, vars_] :=  (* {{{ *)
 Module[{cl, center, Aq, Bq, Cq, Dq, Eq, Fq, cl2, Am},
  cl = CoefficientList[poly, vars];
  {Aq = cl[[3, 1]], Bq = cl[[2, 2]]/2, Cq = cl[[1, 3]], 
   Dq = cl[[2, 1]]/2, Eq = cl[[1, 2]]/2, Fq = cl[[1, 1]]}; 
  center = -Inverse[{{Aq, Bq}, {Bq, Cq}}].{Dq, Eq};
  cl2 = CoefficientList[poly /. {vars[[1]] -> vars[[1]] + center[[1]], vars[[2]] -> vars[[2]] + center[[2]]}, vars];
  Am = {{cl2[[3, 1]], cl2[[2, 2]]/2}, {cl2[[2, 2]]/2, cl2[[1, 3]]}};
  {center, ArcTan @@ (Eigenvectors[Am][[1]]), 
   1/Sqrt[-Eigenvalues[Am]/cl2[[1, 1]]]}
  ]/; PolynomialQ[poly, vars] (* }}} *)
(*  *)
DistanceBetweenSetsOfPoints[p1_List, p2_List] /; 
  If[Length[p1] == Length[p2], True, 
   Message[DistanceBetweenSetsOfPoints::nnarg, Length[p1], 
    Length[p2]]; False] := 
 Module[{p2tmp = p2, n = Length[p2], OrderedList = {}, k},
  Do[
   k = Nearest[p2tmp[[;; n + 1 - i]] -> Automatic, p1[[i]]];
   OrderedList = Append[OrderedList, p2tmp[[k]][[1]]];
   p2tmp = Drop[p2tmp, k];, {i, n}];
  Total[EuclideanDistance @@@ Transpose[{p1, OrderedList}]]]
DistanceBetweenSetsOfPoints::nnarg = 
  "The lengths of the lists must be equal, but they are  `1` and \
`2`.";

HistogramListPoints[data_, Options___] :=Transpose[{Drop[(#[[1]] + RotateLeft[#[[1]]])/
      2, -1], #[[2]]} &[HistogramList[data, Options]]]


HistogramPointsForLine[data_, Options___] := 
 Module[{hrs = HistogramList[data, Options]},
  Transpose[{Flatten[Transpose[{hrs[[1]], hrs[[1]]}]], 
    Flatten[{0, Transpose[{hrs[[2]], hrs[[2]]}], 0}]}]]


RandomUnitVector[n_] := Module[{v},
  v = RandomReal[NormalDistribution[0, 1], n];
  v/Norm[v]]
RandomUnitVector[] := RandomUnitVector[3]


(*
Se comenta porque ya Mathematica trae una funcion que hace eso
(* From http://mathworld.wolfram.com/BlockDiagonalMatrix.html*)
BlockDiagonalMatrix[b : {__?MatrixQ}] := 
 Module[{r, c, n = Length[b], i, j}, {r, c} = 
   Transpose[Dimensions /@ b];
  ArrayFlatten[
   Table[If[i == j, b[[i]], ConstantArray[0, {r[[i]], c[[j]]}]], {i, 
     n}, {j, n}]]]
*)

NumberList[lista_]:=Flatten[Evaluate[#], 1] & /@ Transpose[{Range[Length[lista]], lista}]

OffSpellErrors[]:={Off[General::spell],Off[General::spell1]}
OnSpellErrors[]:={On[General::spell],On[General::spell1]}

ReadListUncomment[file_, Options___] := 
 ReadList[ StringToStream[
   StringJoin[ StringInsert[#, "\n", -1] & /@  Select[ReadList[file, String], StringFreeQ[#, "#"] &]]], Options]


ColumnAddKeepFirst[MultiList_] :=  MapThread[Prepend, {(Plus @@ MultiList)[[All, 2 ;;]], MultiList[[1, All, 1]]}]
ColumnAddKeepFirst[FirstList_, SecondList_] := ColumnAddKeepFirst[{FirstList, SecondList}]

(*  Legends and symbols {{{ *)
InsetWithSymbols[LowerLeft_List,BoxSize_List,RealtiveCoordinateLowerSymbol_, SepSymbols_,SymbolList_,TextList_, TextSpacing_]:=
	Module[{i},
	{Table[ SymbolList[[i]][ LowerLeft+RealtiveCoordinateLowerSymbol+{0,(i-1) SepSymbols}],{i, Length[SymbolList]}],
      Graphics[ Table[Text[TextList[[i]], LowerLeft+ RealtiveCoordinateLowerSymbol+{0,(i-1) SepSymbols}+{TextSpacing, 0},{-1,0}]
                               ,{i,Length[TextList]}]],
      Graphics[{Line[{LowerLeft,LowerLeft+{0,BoxSize[[2]]},LowerLeft+BoxSize, LowerLeft+{BoxSize[[1]],0},LowerLeft}]}]}
]


MyTriangle[{x_, y_}, Color1_, Proportion_, delta_, th_] := 
    Graphics[{Color1,Polygon[{Scaled[delta {-1/2, -Proportion/3}, {x, y}], 
       Scaled[delta {0, 2 Proportion/3}, {x, y}],Scaled[delta {1/2, -Proportion/3}, {x, y}]}], Thickness[th], 
       GrayLevel[0], Line[{Scaled[delta {-1/2, -Proportion/3}, {x, y}], 
       Scaled[delta {0, 2 Proportion/3}, {x, y}], Scaled[delta {1/2, -Proportion/3}, {x, y}], 
       Scaled[delta {-1/2, -Proportion/3}, {x, y}]}]}];

MyInvertedTriangle[{x_, y_}, Color1_, Proportion_, delta_, th_] := 
    Graphics[{Color1,Polygon[{Scaled[delta {-1/2, Proportion/3}, {x, y}], 
       Scaled[delta {0, -2 Proportion/3}, {x, y}],Scaled[delta {1/2, Proportion/3}, {x, y}]}], Thickness[th], 
       GrayLevel[0], Line[{Scaled[delta {-1/2, Proportion/3}, {x, y}], 
       Scaled[delta {0, -2 Proportion/3}, {x, y}], Scaled[delta {1/2, Proportion/3}, {x, y}], 
       Scaled[delta {-1/2, Proportion/3}, {x, y}]}]}];

MySquare[{x_, y_}, Color1_, Proportion_, delta_, th_] := 
    Graphics[{Color1,  Rectangle[Scaled[delta{-1/2, -Proportion/2}, {x, y}], 
      Scaled[delta{1/2, Proportion/2}, {x, y}]], Thickness[th],GrayLevel[0], 
      Line[{Scaled[delta{-1/2, -Proportion/2}, {x, y}],Scaled[delta{-1/2, Proportion/2}, {x, y}], 
      Scaled[delta{1/2, Proportion/2}, {x, y}], Scaled[delta{1/2, -Proportion/2}, {x, y}], 
      Scaled[delta{-1/2, -Proportion/2}, {x, y}]}]}];

MyRhombous[{x_, y_}, Color1_, Proportion_, delta_, th_] := 
    Graphics[{Color1, 
        Polygon[{Scaled[delta{0, -Proportion/2}, {x, y}], 
            Scaled[delta{1/2, 0}, {x, y}], 
            Scaled[delta{0, Proportion/2}, {x, y}], 
            Scaled[delta{-1/2, 0}, {x, y}]}], Thickness[th], GrayLevel[0], 
        Line[{Scaled[delta{0, -Proportion/2}, {x, y}], 
            Scaled[delta{1/2, 0}, {x, y}], 
            Scaled[delta{0, Proportion/2}, {x, y}], 
            Scaled[delta{-1/2, 0}, {x, y}],  
            Scaled[delta{0, -Proportion/2}, {x, y}]}]}];

My4PointStar1[{x_, y_}, Color1_, Proportion_, delta_, th_] := 
  Module[{PointSet, PointSetLine, theta, alpha},
    alpha = .2;
    PointSet = {{1, 0}, alpha{1, 1}, {0, 1}, alpha{-1, 1}, {-1, 0}, 
        alpha{-1, -1}, {0, -1}, alpha{1, -1}};
    PointSetLine = Flatten[{PointSet, {PointSet[[1]]}}, 1];
    Graphics[{Color1, 
        Polygon[Scaled[delta {#[[1]], Proportion #[[2]]}, {x, y}] & /@ 
            PointSet], Thickness[th], GrayLevel[0], 
        Line[Scaled[delta {#[[1]], Proportion #[[2]]}, {x, y}] & /@ 
            PointSetLine]}]]

My4PointStar2[{x_, y_}, Color1_, Proportion_, delta_, th_] := 
  Module[{PointSet, PointSetLine, theta, alpha},
    alpha = .3;
    PointSet = {alpha{1, 0}, {1, 1}, alpha{0, 1}, {-1, 1}, 
        alpha{-1, 0}, {-1, -1}, alpha{0, -1}, {1, -1}};
    PointSetLine = Flatten[{PointSet, {PointSet[[1]]}}, 1];
    Graphics[{Color1, 
        Polygon[Scaled[delta {#[[1]], Proportion #[[2]]}, {x, y}] & /@ 
            PointSet], Thickness[th], GrayLevel[0], 
        Line[Scaled[delta {#[[1]], Proportion #[[2]]}, {x, y}] & /@ 
            PointSetLine]}]]


My5PointStar[{x_, y_}, Color1_, Proportion_, delta_, th_] := Module[{PointSet, PointSetLine, theta,theta2}, 
      PointSet = Flatten[Table[
      {{Cos[theta + Pi/2], Sin[theta + Pi/2]}, 
              1/2 (3 - Sqrt[5]) {Cos[theta + Pi/5 + Pi/2], 
                  Sin[theta + Pi/5 + Pi/2]}}, {theta, 0, 2  Pi - 2  Pi/5, 2  Pi/5}], 1];
     PointSetLine =  Flatten[{PointSet, {PointSet[[1]]}}, 1];
     Graphics[{Color1, Polygon[Scaled[delta {#[[1]], Proportion #[[2]]}, {x, y}] & /@ PointSet], 
        Thickness[th], GrayLevel[0], 
          Line[Scaled[delta {#[[1]], Proportion #[[2]]}, {x, y}] &/@ PointSetLine]}]]

MyInverted5PointStar[{x_, y_}, Color1_, Proportion_, delta_, th_] := Module[{PointSet, PointSetLine, theta,theta2}, 
      PointSet = Flatten[Table[
      theta=theta2+Pi;
      {{Cos[theta + Pi/2], Sin[theta + Pi/2]}, 
              1/2 (3 - Sqrt[5]) {Cos[theta + Pi/5 + Pi/2], 
                  Sin[theta + Pi/5 + Pi/2]}}, {theta2, 0, 2  Pi - 2  Pi/5, 2  Pi/5}], 1];
     PointSetLine =  Flatten[{PointSet, {PointSet[[1]]}}, 1];
     Graphics[{Color1, Polygon[Scaled[delta {#[[1]], Proportion #[[2]]}, {x, y}] & /@ PointSet], 
        Thickness[th], GrayLevel[0], 
          Line[Scaled[delta {#[[1]], Proportion #[[2]]}, {x, y}] &/@ PointSetLine]}]]

MyCircle[{x_, y_}, Color1_, Proportion_, delta_, th_] := 
  Graphics[{Color1, Disk[{x, y}, Scaled[delta{1/2, Proportion/2}]], 
      Thickness[th], GrayLevel[0], Circle[{x, y}, Scaled[delta{1/2, Proportion/2}]]}]

MyEllipse1[{x_, y_}, Color1_, Proportion_, delta_, th_] := 
  Graphics[{Color1, Disk[{x, y}, Scaled[delta{.5 1/2, Proportion/2 8/5}]], 
      Thickness[th], GrayLevel[0], Circle[{x, y}, Scaled[delta{.5 1/2, Proportion/2 8/5}]]}]
MyEllipse2[{x_, y_}, Color1_, Proportion_, delta_, th_] := 
  Graphics[{Color1, Disk[{x, y}, Scaled[delta{.8, .5 Proportion/2}]], 
      Thickness[th], GrayLevel[0], Circle[{x, y}, Scaled[delta{.8, .5 Proportion/2}]]}]

My4PointStar3[{x_, y_}, Color1_, Proportion_, delta_, th_] := 
 Module[{PointSet, PointSetLine, theta, alpha}, alpha = .3;
  PointSet = {{0, 0}, {Cos[\[Pi]/4 - alpha], 
     Sin[\[Pi]/4 - alpha]}, {Cos[\[Pi]/4 + alpha], 
     Sin[\[Pi]/4 + alpha]}, {0, 0}, {Cos[3 \[Pi]/4 - alpha], 
     Sin[3 \[Pi]/4 - alpha]}, {Cos[3 \[Pi]/4 + alpha], 
     Sin[3 \[Pi]/4 + alpha]}, {0, 0}, {Cos[5 \[Pi]/4 - alpha], 
     Sin[5 \[Pi]/4 - alpha]}, {Cos[5 \[Pi]/4 + alpha], 
     Sin[5 \[Pi]/4 + alpha]}, {0, 0}, {Cos[7 \[Pi]/4 - alpha], 
     Sin[7 \[Pi]/4 - alpha]}, {Cos[7 \[Pi]/4 + alpha], 
     Sin[7 \[Pi]/4 + alpha]}};
  PointSetLine = Flatten[{PointSet, {PointSet[[1]]}}, 1];
  Graphics[{Color1, 
    Polygon[Scaled[delta {#[[1]], Proportion #[[2]]}, {x, y}] & /@ 
      PointSet], Thickness[th], GrayLevel[0], 
    Line[Scaled[delta {#[[1]], Proportion #[[2]]}, {x, y}] & /@ 
      PointSetLine]}]]
MyLegend[TheStyle_List, Heigth_, Xpos_, Xlength_, TheText_, XSepText_,  OptionsPattern[]] := 
  {Text[TheText, Scaled[{Xpos + Xlength + XSepText, Heigth}],OptionValue[Alignment]], 
  Join[TheStyle, {Line[{Scaled[{Xpos, Heigth}], 
      Scaled[{Xpos + Xlength, Heigth}]}]}]}

LegendBox[TheLegends_, TheStyles_, UpperHeight_, Xpos_, Xlength_, 
  XSepText_, \[CapitalDelta]Height_,  OptionsPattern[]] := 
 Module[{i}, 
  Table[MyLegend[TheStyles[[i]], 
    UpperHeight - (i - 1) \[CapitalDelta]Height, Xpos, Xlength, 
    TheLegends[[i]], XSepText, Alignment -> OptionValue[Alignment]], {i, Length[TheLegends]}]]
Options[LegendBox] = {Alignment->{0,0}};
Options[MyLegend] = {Alignment->{0,0}};


MySymbol[Coordinate_,  OptionsPattern[]] :=
 {MyTriangle,MySquare,MyRhombous,MyInvertedTriangle,MyCircle,My5PointStar,
 	My4PointStar1,My4PointStar2,MyInverted5PointStar, My4PointStar3, MyEllipse1,MyEllipse2}[[OptionValue[SymbolNumber]]][Coordinate, 
  OptionValue[Color], OptionValue[Proportion], OptionValue[delta], OptionValue[ThicknessBorder]]
Options[MySymbol] = {SymbolNumber -> 1, Color -> Hue[0],
	Proportion -> GoldenRatio, delta -> 0.02,  ThicknessBorder -> 0.001};
(*  }}}  *)

MyAverage[x_] := Plus @@ x/Length[x]

Seed[] := Floor[Random[] 1000000000]

Instruction[jodas_List, Executable_String] := Module[{tmpins},
      tmpins = "printf \""; 
      Do[Do[tmpins = tmpins <> ToString[jodas[[i, j]]] <> " ";, {j, 
            Length[ jodas[[i]] ]}]; 
        tmpins = tmpins <> "\\n";, {i, Length[jodas]}];
      tmpins <> "\" | " <> Executable];

Norma[x_List] := Sqrt[Plus @@( (Abs[x])^2)]


ColorCoding[NumberOfNumbers_Integer,NumberOfColors_Integer]:=
  Module[{n1,n2},n1=2 Pi/NumberOfNumbers;
    n2=2 Pi/NumberOfColors;
    Show[{Graphics[({Hue[#1/(2 Pi-n1)],
                  Text[ToString[N[#1/(2. Pi),2]],1.2 {Cos[#1],Sin[#1]}]}&)/@
            Range[0,2 Pi-n1,n1]],
        Graphics[({Hue[#1/(2 Pi-n2)],Disk[{0,0},1,{#1,#1+n2}]}&)/@
            Range[0,2 Pi-n2,n2]]},DisplayFunction\[Rule]Identity,
      AspectRatio\[Rule]Automatic]]

    
End[];
EndPackage[];
