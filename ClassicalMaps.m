(* {{{ *) BeginPackage["ClassicalMaps`"]
(* {{{ Maps *)
sawTooth::usage="sawTooth[{q_, p_}, K_] gives the sawtooth as in DOI:	10.1016/S0375-9601(97)00455-6"
harper::usage="Harper map, as in our article with Diego Wisniacki, usage: harper[{q_, p_}, {k_, kp_}], harper[{q_, p_}, k_]"
standardMap::usage="Standard Map, [0,1]x[0,1] standardMap[{q_, p_}, K_]"
standardMapNonPeriodic::usage="Standard Map, without periodic conditions, [0,1]x[0,1] standardMapNonPeriodic[{q_, p_}, K_]"
(* }}} *)
Begin["`Private`"] 
(* {{{ Maps *)
(* {{{ *) sawTooth[{q_, p_}, K_] := Module[{qp, pp},
  pp = p + K (Mod[q, 1] - .5);
  qp = q + pp;
  {qp, pp}] 
(* }}} *)
(* {{{ *) standardMap[{q_, p_}, K_] :=
	{Mod[#[[1]],1],#[[2]]}&[standardMapNonPeriodic[{q, p}, K]]
(* }}} *)
(* {{{ *) standardMapNonPeriodic[{q_, p_}, K_] := 
	Module[{qp, pp}, 
		pp = p + K (Sin[2 \[Pi] q])/(2 \[Pi]);
  		qp = q + pp;
 	 {qp, pp}]



(* }}} *)
(* {{{ *) harper[{q_, p_}, {k_, kp_}] := Module[{qp, pp},
  pp = p - k Sin[2. \[Pi] q];
  qp = q + kp Sin[2. \[Pi] pp];
  {qp, pp}]
(* }}} *) 
(* {{{ *) harper[{q_, p_}, k_] := harper[{q, p}, {k, k}]
(* }}} *)
(* {{{ *) 
(* }}} *)
(* }}} *)
End[] 
EndPackage[] (* }}} *)
