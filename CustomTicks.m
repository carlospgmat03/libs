(* :Title: CustomTicks *)
(* :Context: CustomTicks` *)
(* :Author: Mark A. Caprio, Center for Theoretical Physics, 
  Yale University *)
(* :Summary: Custom tick mark generation for linear, log, 
  and general nonlinear axes. *)
(* :Copyright: Copyright 2005, Mark A. Caprio *)
(* :Package Version: 1.2 *)
(* :Mathematica Version: 4.0 *)
(* :History:
      MCAxes package, January 10, 2003.
     MCAxes and then MCTicks packages distributed as part of LevelScheme, 
  2004.
     V1.0, March 11, 2005.  MathSource No. 5599.
     V1.1, March 18, 2005.  Documentation update.
      V1.2, September 17, 2005.  Simplified LogTicks syntax.
  *)



BeginPackage["CustomTicks`"];

Unprotect[Evaluate[$Context<>"*"]];





LinTicks::usage="LinTicks[x1,x2,spacing,subdivisions] or LinTicks[x1,x2] or LinTicks[majorlist,minorlist].";\

LogTicks::usage="LogTicks[power1,power2,subdivisions] or LogTicks[base,power1,power2,subdivisions].";\


TickPreTransformation::usage="Option for LinTicks.  Mapping from given coordinate value to value used for range tests labeling.";\

TickPostTransformation::usage="Option for LinTicks.  Mapping from tick values to actual coordinates used for positioning tick mark, applied after all range tests and labeling.";\

TickRange::usage="Option for LinTicks.";

ShowTickLabels::usage="Option for LinTicks.";
TickLabelRange::usage="Option for LinTicks.";
ShowFirst::usage="Option for LinTicks.";
ShowLast::usage="Option for LinTicks.";
ShowMinorTicks::usage="Option for LinTicks.";
TickLabelStart::usage="Option for LinTicks.";
TickLabelStep::usage="Option for LinTicks.";
TickLabelFunction::usage="Option for LinTicks.";
DecimalDigits::usage="Option for LinTicks.";

MajorTickLength::usage="Option for LinTicks.";
MinorTickLength::usage="Option for LinTicks.";
MajorTickStyle::usage="Option for LinTicks.";
MinorTickStyle::usage="Option for LinTicks.";
MinorTickIndexRange::usage="Option for LinTicks.";
MinorTickIndexTransformation::usage="Option for LinTicks.";

FractionDigits::usage="FractionDigits[x] returns the number of digits to the right of the point in the decimal representation of x.  It will return large values, determined by Precision, for some numbers, e.g., non-terminating rationals.";\

FractionDigitsBase::usage="Option for FractionDigits.";

LimitTickRange::usage="LimitTickRange[{x1,x2},ticks] selects those ticks with coordinates approximately in the range x1...x2.  Ticks must be specified in full form, or at least as a list rather than a bare number.";\

TransformTicks::usage=
    "StripTickLabels[positionfcn,lengthfcn,ticks] transforms the positions and lengths of all tick marks in a list.  Tick marks must be specified in full form, or at least with an explicit pair of in and out lengths.";\

StripTickLabels::usage="StripTickLabels[ticks] removes any text labels from ticks.  Ticks must be specified in full form, or at least as a list rather than a bare number.";\



AugmentTicks::usage="AugmentTicks[labelfunction,lengthlist,stylelist,ticks] augments any ticks in ticklist to full form.";\

AugmentAxisTickOptions::usage="AugmentAxisTickOptions[numaxes,list] replaces any None entries with null lists and appends additional null lists as needed to make numaxes entries.  AugmentAxisTickOptions[numaxes,None] returns a list of null lists.  Note that this differs from the behavior of the Mathematica plotting functions with FrameTicks, for which an unspecified upper or right axis duplicates the given lower or left axis.";\

TickQ::usage="TickQ[x] returns True if x is a valid tick specification." ;
TickListQ::usage="TickListQ[x] returns True if x is a list of valid tick specifications.  This is not a general test for a valid axis tick option specification, as None, Automatic, and a function can also be valid axis tick option specifications.";\
 



LogTicks::oldsyntax="The number of minor subdivisions no longer needs to be specified for LogTicks (see CustomTicks manual for new syntax).";
LogTicks::minorsubdivs="Number of minor subdivisions `1` specified for LogTicks is not 1 or \[LeftCeiling]base\[RightCeiling]-1 (i.e., \[LeftCeiling]base\[RightCeiling]-2 tick marks) and so is being ignored.";\

AugmentTicks::automatic="Tick list must be specified explicitly.";
AugmentAxisTickOptions::numaxes="Tick lists specified for more than `1` axes.";



Begin["`Private`"];



(* range testing and manipulation utilities, from package MCGraphics *)

InRange[{x1_,x2_},x_]:=(x1\[LessEqual]x)&&(x\[LessEqual]x2);

InRange[{x1_,x2_},x_]:=(x1\[LessEqual]x)&&(x\[LessEqual]x2);
InRangeProper[{x1_,x2_},x_]:=(x1<x)&&(x<x2);

InRegion[{{x1_,x2_},{y1_,y2_}},{x_,y_}]:=
    InRange[{x1,x2},x]&&InRange[{y1,y2},y];
InRegion[{x1_,y1_},{x2_,y2_},{x_,y_}]:=InRange[{x1,x2},x]&&InRange[{y1,y2},y];

ExtendRange[PRange:{x1_,x2_},PFrac:{fx1_,fx2_}]:=
    PRange+PFrac*{-1,+1}*-Subtract@@PRange;
ExtendRegion[PRange:{{x1_,x2_},{y1_,y2_}},PFrac:{{fx1_,fx2_},{fy1_,fy2_}}]:=
    PRange+PFrac*{{-1,+1},{-1,+1}}*-Subtract@@@PRange;

(* approximate equality testing utility, from package MCArithmetic *)

Options[ApproxEqual]={Chop\[Rule]1*^-10};
ApproxEqual[x_?NumericQ,y_?NumericQ,Opts___?OptionQ]:=Module[
      {FullOpts=Flatten[{Opts,Options[ApproxEqual]}]},
      Chop[x-y,Chop/.FullOpts]\[Equal]0
      ];
ApproxEqual[x_?NumericQ,_DirectedInfinity]:=False;
ApproxEqual[_DirectedInfinity,x_?NumericQ]:=False;
ApproxInRange[{x1_,x2_},x_,Opts___?OptionQ]:=
    ApproxEqual[x,x1,Opts]||((x1\[LessEqual]x)&&(x\[LessEqual]x2))||
      ApproxEqual[x,x2,Opts];
ApproxIntegerQ[x_?NumericQ,Opts___?OptionQ]:=ApproxEqual[Round[x],x,Opts];







Options[FractionDigits]={FractionDigitsBase\[Rule]10,Limit\[Rule]Infinity,
      Chop\[Rule]1*^-8};

FractionDigits[x_?NumericQ,Opts___?OptionQ]:=Module[
      {
        FullOpts=Flatten[{Opts,Options[FractionDigits]}],
        Value,NumToRight,OptFractionDigitsBase,OptLimit
        },
      Value=N[x];
      OptFractionDigitsBase=FractionDigitsBase/.FullOpts;
      OptLimit=Limit/.FullOpts;
      NumToRight=0;
      While[
        !ApproxIntegerQ[Value,Chop\[Rule](Chop/.FullOpts)]&&(NumToRight<
              OptLimit),
        Value*=OptFractionDigitsBase;
        NumToRight++
        ];
      NumToRight
      ];





NumberSignsOption[x_]:=Switch[
      Sign[x],
      0|(+1),NumberSigns\[Rule]{"",""},
      -1,NumberSigns\[Rule]{"-",""}
      ];

FixedTickFunction[ValueList_List,DecimalDigits_]:=Module[
      {PaddingDecimals,
        PaddingDigitsLeft,PaddingDigitsRight
        },
      PaddingDigitsLeft=If[
          ValueList==={},
          0,
          Length[IntegerDigits[IntegerPart[Max[Abs[ValueList]]]]]
          ];
      PaddingDigitsRight=If[
          DecimalDigits===Automatic,
          If[
            ValueList==={},
            0,
            Max[FractionDigits/@ValueList]
            ],
          DecimalDigits
          ];
      If[PaddingDigitsRight\[Equal]0,
        PaddedForm[IntegerPart[#],PaddingDigitsLeft,NumberSignsOption[#]]&,
        PaddedForm[#,{PaddingDigitsLeft+PaddingDigitsRight,
              PaddingDigitsRight},NumberSignsOption[#]]&
        ]
      ];





Options[LinTicks]={TickPreTransformation\[Rule]Identity,
      TickPostTransformation\[Rule]Identity,ShowFirst\[Rule]True,
      ShowLast\[Rule]True,ShowTickLabels\[Rule]True,ShowMinorTicks\[Rule]True,
      TickLabelStart\[Rule]0,TickLabelStep\[Rule]1,
      TickRange\[Rule]{-Infinity,Infinity},
      TickLabelRange\[Rule]{-Infinity,Infinity},
      TickLabelFunction\[Rule]Automatic,DecimalDigits\[Rule]Automatic,
      MajorTickLength\[Rule]{0.010,0},MinorTickLength\[Rule]{0.005,0},
      MajorTickStyle\[Rule]{},MinorTickStyle\[Rule]{},
      MinorTickIndexRange\[Rule]{1,Infinity},
      MinorTickIndexTransformation\[Rule]Identity};

LinTicks[RawMajorCoordList_List,RawMinorCoordList_List,Opts___]:=Module[
      {
        FullOpts= Flatten[{Opts,Options[LinTicks]}],
        MajorCoordList,
        LabeledCoordList,
        MinorCoordList,
        MajorTickList,
        MinorTickList,
        UsedTickLabelFunction,
        DefaultTickLabelFunction,
        TickValue,
        TickPosition,
        TickLabel,
        TickLength,
        TickStyle,
        i
        },
      
      (* make major ticks *)
      MajorCoordList=
        Select[(TickPreTransformation/.FullOpts)/@RawMajorCoordList,
          ApproxInRange[(TickRange/.FullOpts),#]&];
      LabeledCoordList=Flatten[Table[
            TickValue=MajorCoordList[[i]];
            If[
              (ShowTickLabels/.FullOpts)
                &&ApproxInRange[(TickLabelRange/.FullOpts),TickValue]
                &&(Mod[
                      i-1,(TickLabelStep/.FullOpts)]\[Equal](TickLabelStart/.\
FullOpts))
                &&((i\[NotEqual]1)||(ShowFirst/.FullOpts))
                &&((i\[NotEqual]
                        Length[MajorCoordList])||(ShowLast/.FullOpts)),
              TickValue,
              {}
              ],
            {i,1,Length[MajorCoordList]}
            ]];
      DefaultTickLabelFunction=
        FixedTickFunction[LabeledCoordList,DecimalDigits/.FullOpts];
      UsedTickLabelFunction=Switch[
          (TickLabelFunction/.FullOpts),
          Automatic,(#2&),
          _,(TickLabelFunction/.FullOpts)
          ];
      TickLength=(MajorTickLength/.FullOpts);
      TickStyle=(MajorTickStyle/.FullOpts);
      MajorTickList=Table[
          
          (* calculate tick value *)
          TickValue=MajorCoordList[[i]];
          
          (* calculate coordinate for drawing tick *)
          TickPosition=(TickPostTransformation/.FullOpts)[TickValue];
          
          (* construct label, 
            or null string if it should be suppressed -- if: major tick, 
            
            in designated modular cycle if only a cycle of major ticks are to \
be labeled, tick is in TickLabelRange, 
            and not explicitly suppressed as first or last label; 
            will only then be used if tick is also in TickRange *)
          TickLabel=If[
              (ShowTickLabels/.FullOpts)
                &&ApproxInRange[(TickLabelRange/.FullOpts),TickValue]
                &&(Mod[
                      i-1,(TickLabelStep/.FullOpts)]\[Equal](TickLabelStart/.\
FullOpts))
                &&((i\[NotEqual]1)||(ShowFirst/.FullOpts))
                &&((i\[NotEqual]Length[
                          MajorCoordList])||(ShowLast/.FullOpts)),
              
              UsedTickLabelFunction[TickValue,
                DefaultTickLabelFunction[TickValue]],
              ""
              ];
          
          (* make tick *)
          {TickPosition,TickLabel,TickLength,TickStyle},
          
          {i,1,Length[MajorCoordList]}
          ];
      
      (* make minor ticks *)
      MinorCoordList=
        Select[(TickPreTransformation/.FullOpts)/@RawMinorCoordList,
          ApproxInRange[(TickRange/.FullOpts),#]&];
      TickLength=(MinorTickLength/.FullOpts);
      TickStyle=(MinorTickStyle/.FullOpts);
      MinorTickList=If[(ShowMinorTicks/.FullOpts),
          Table[
            
            (* calculate tick value *)
            TickValue=MinorCoordList[[i]];
            
            (* calculate coordinate for drawing tick *)
            TickPosition=(TickPostTransformation/.FullOpts)[TickValue];
            
            (* make tick *)
            {TickPosition,"",TickLength,TickStyle},
            
            {i,1,Length[MinorCoordList]}
            ],
          {}
          ];
      
      (* combine tick lists*)
      Join[MajorTickList,MinorTickList]
      
      ];

LinTicks[x1_?NumericQ,x2_?NumericQ,Opts___?OptionQ]:=Module[
    {UsedRange,DummyGraphics,TickList,MajorCoordList,MinorCoordList,x},
    
    (* extend any round-number range by a tiny amount *)
    (* this seems to make Mathematica 4.1 give a much cleaner, 
      sparser set of ticks *)
    UsedRange=If[
        And@@ApproxIntegerQ/@{x1,x2},
        ExtendRange[{x1,x2},{1*^-5,1*^-5}],
        {x1,x2}
        ]; 
    
    (* extract raw tick coordinates from Mathematica *)
    DummyGraphics=
      Show[Graphics[{}],PlotRange\[Rule]{UsedRange,Automatic},
        DisplayFunction\[Rule]Identity];
    TickList=First[Ticks/.AbsoluteOptions[DummyGraphics,Ticks]];
    MajorCoordList=Cases[TickList,{x_,_Real,___}\[RuleDelayed]x];
    MinorCoordList=Cases[TickList,{x_,"",___}\[RuleDelayed]x];
    
    (* generate formatted tick mark specifications *)
    LinTicks[MajorCoordList,MinorCoordList,Opts]
    ]

LinTicks[x1_?NumericQ,x2_?NumericQ,Spacing_?NumericQ,MinorSubdivs:_Integer,
      Opts___]:=Module[
      {
        FullOpts= Flatten[{Opts,Options[LinTicks]}],
        MaxMajorIndex,
        MajorCoordList,MinorCoordList
        },
      
      (* preliminary calculations  *)
      MaxMajorIndex=Round[(x2-x1)/Spacing];
      
      (* construct table of ticks --indexed by MajorIndex=0,1,...,
        MaxMajorTick and MinorIndex=0,...,MinorSubdivs-1, 
        where MinorIndex=0 gives the major tick, 
        except no minor ticks after last major tick *)
      MajorCoordList=Flatten[Table[
            N[x1+MajorIndex*Spacing+MinorIndex*Spacing/MinorSubdivs],
            {MajorIndex,0,MaxMajorIndex},{MinorIndex,0,0}
            ]];
      MinorCoordList=Flatten[Table[
            If[
              InRange[MinorTickIndexRange/.FullOpts,MinorIndex],
              
              N[x1+MajorIndex*
                    Spacing+((MinorTickIndexTransformation/.FullOpts)@
                        MinorIndex)*Spacing/MinorSubdivs],
              {}
              ],
            {MajorIndex,0,MaxMajorIndex},{MinorIndex,1,MinorSubdivs-1}
            ]];
      
      (* there are usually ticks to be suppressed at the upper end, 
        since the major tick index rounds up to the next major tick (for \
safety in borderline cases where truncation might fail), 
        and the loop minor tick index iterates for a full series of minor \
ticks even after the last major tick *)
      MajorCoordList=Select[MajorCoordList,ApproxInRange[{x1,x2},#]&];
      MinorCoordList=Select[MinorCoordList,ApproxInRange[{x1,x2},#]&];
      
      LinTicks[MajorCoordList,MinorCoordList,Opts]
      
      ];





LogTicks[Base:(_?NumericQ):10,p1_Integer,p2_Integer,Opts___?OptionQ]:=Module[
      {BaseSymbol,MinorSubdivs},
      
      BaseSymbol=If[Base===E,StyleForm["e",FontFamily->"Italic"],Base];
      MinorSubdivs=Ceiling[Base]-1; (* one more than minor ticks *)
      MinorSubdivs=Max[MinorSubdivs,1]; (* 
        prevent underflow from bases less than 2 *)
      LinTicks[p1,p2,1,MinorSubdivs,
        TickLabelFunction\[Rule](DisplayForm[
                SuperscriptBox[BaseSymbol,IntegerPart[#]]]&),
        MinorTickIndexTransformation\[Rule](Log[Base,#+1]*MinorSubdivs&),
        Opts
        ]
      
      ];

(* syntax traps for old syntax -- but will not catch usual situation in which \
base was unspecified but subdivs was *)
LogTicks[Base_?NumericQ,p1_Integer,p2_Integer,MinorSubdivs_Integer,
        Opts___?OptionQ]/;(MinorSubdivs\[Equal]
          Max[Ceiling[Base]-1,1]):=(Message[LogTicks::oldsyntax];
      LogTicks[Base,p1,p2,ShowMinorTicks\[Rule]True,Opts]);
LogTicks[Base_?NumericQ,p1_Integer,p2_Integer,MinorSubdivs_Integer,
        Opts___?OptionQ]/;(MinorSubdivs==1):=(Message[LogTicks::oldsyntax];
      LogTicks[Base,p1,p2,ShowMinorTicks\[Rule]False,Opts]);
LogTicks[Base_?NumericQ,p1_Integer,p2_Integer,MinorSubdivs_Integer,
        Opts___?OptionQ]/;((MinorSubdivs!=
              Max[Ceiling[Base]-1,1])&&(MinorSubdivs!=1)):=(Message[
        LogTicks::oldsyntax];Message[LogTicks::minorsubdivs,MinorSubdivs];
      LogTicks[Base,p1,p2,ShowMinorTicks\[Rule]True,Opts]);





AugmentTick[LabelFunction_,DefaultLength_List,DefaultStyle_List,x_?NumericQ]:=
    AugmentTick[LabelFunction,DefaultLength,DefaultStyle,{x}];



AugmentTick[LabelFunction_,DefaultLength_List,
      DefaultStyle_List,{x_?NumericQ}]:=
    AugmentTick[LabelFunction,DefaultLength,DefaultStyle,{x,LabelFunction@x}];



AugmentTick[LabelFunction_,DefaultLength_List,
      DefaultStyle_List,{x_?NumericQ,LabelText_}]:=
    AugmentTick[LabelFunction,DefaultLength,
      DefaultStyle,{x,LabelText,DefaultLength}];



AugmentTick[LabelFunction_,DefaultLength_List,
      DefaultStyle_List,{x_?NumericQ,LabelText_,PLen_?NumericQ,RestSeq___}]:=
    AugmentTick[LabelFunction,DefaultLength,
      DefaultStyle,{x,LabelText,{PLen,0},RestSeq}];



AugmentTick[LabelFunction_,DefaultLength_List,
      DefaultStyle_List,{x_?NumericQ,LabelText_,LengthList_List}]:=
    AugmentTick[LabelFunction,DefaultLength,
      DefaultStyle,{x,LabelText,LengthList,DefaultStyle}];



AugmentTick[LabelFunction_,DefaultLength_List,DefaultStyle_List,
      TheTick:{x_?NumericQ,LabelText_,LengthList_List,Style_List}]:=TheTick;



AugmentTicks[DefaultLength_List,DefaultStyle_List,TickList_List]:=
    AugmentTick[""&,DefaultLength,DefaultStyle,#]&/@TickList;





AugmentAxisTickOptions[NumAxes_Integer,TickLists:None]:=Table[{},{NumAxes}];
AugmentAxisTickOptions[NumAxes_Integer,TickLists_List]:=Module[
      {},
      If[
        NumAxes<Length[TickLists],
        Message[AugmentAxisTickOptions::numaxes,NumAxes]
        ];
      Join[
        Replace[TickLists,{None\[Rule]{}},{1}],
        Table[{},{NumAxes-Length[TickLists]}]
        ]
      ];





TickPattern=(_?NumericQ)|{_?NumericQ}|{_?NumericQ,_}|{_?
          NumericQ,_,(_?NumericQ)|{_?NumericQ,_?NumericQ}}|{_?
          NumericQ,_,(_?NumericQ)|{_?NumericQ,_?NumericQ},_List};

TickQ[x_]:=MatchQ[x,TickPattern];
TickListQ[x_]:=MatchQ[x,{}|{TickPattern..}];















LimitTickRange[Range:{x1_,x2_},TickList_List]:=
    Cases[TickList,{x_,___}/;ApproxInRange[{x1,x2},x]];

TransformTicks[PosnTransformation_,LengthTransformation_,TickList_List]:=
    Replace[TickList,{x_,t_,l:{_,_},
          RestSeq___}\[RuleDelayed]{PosnTransformation@x,t,
          LengthTransformation/@l,RestSeq},{1}];

StripTickLabels[TickList_List]:=
    Replace[TickList,{x_,_,RestSeq___}\[RuleDelayed]{x,"",RestSeq},{1}];



End[];



Protect[Evaluate[$Context<>"*"]];
EndPackage[];













































































































































































