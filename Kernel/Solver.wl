(* ::Package:: *)

(* ::Section:: *)
(*Package Headers*)


Begin["Taggar`GERF`Private`"];


(* ::Section:: *)
(*Definitions*)


ClearAll[GERFSolve];
Options[GERFSolve] = {
	"RationalFunction" -> Function[n, 1/(1+E^n)],
	"WaveConstants" -> Automatic,
	"OutputMode" -> "Solutions",
	"BalanceConstant" -> Automatic
};
GERFSolve::usage =
	"GERFSolve[eqn, u[x, t]] solves the given eqn in u[x, t] using GERF expansion technique.
GERFSolve[{eqn1, eqn2, ..., eqnk}, {u1[x, t], ..., uk[x, t]}] solves the given {eqn1, eqn2, ..., eqnj} in u1[x, t], ..., uk[x, t]} using GERF expansion technique.";
GERFSolve::GERFPackageError =
	"`1`";
GERFSolve::InvalidFractionalDerivatives =
	"Multiple fractional orders received for these dimensions: `1`. This is not allowed.";
GERFSolve::ConstantsLengthMismatch =
	"The number of wave constants provided in the \"WaveConstants\" option does not match the number of independent variables in the equation(s).";


GERFSolve[ieqns_List, us : {(_[vars__]) ..}, opts : OptionsPattern[]] :=
	Module[
		{w, A, R, U, eta, state}, (* initialise A,R,U,eta here so each trial solution has same instance *)
		
		state = <|
			"OriginalEquation" -> ieqns,
			"Function" -> AssociationThread[Range @ Length @ ieqns -> Head /@ us],
			"Variables" -> {vars},
			"Length" -> Length @ ieqns,
			"Options" -> <|
				"RationalFunction" -> OptionValue["RationalFunction"],
				"WaveConstants" -> OptionValue["WaveConstants"],
				"OutputMode" -> OptionValue["OutputMode"],
				"BalanceConstant" -> OptionValue["BalanceConstant"]
			|>,
			"ODE" -> None,
			"TrialSolution" -> None,
			"AuxiliaryPolynomial" -> None,
			"AuxiliaryFunction" -> AssociationMap[U, Range @ Length @ ieqns],
			"Eta" -> eta,
			"WaveConstant" -> AssociationMap[w, {vars}],
			"WCH" -> w, (* ad hoc *)
			"BalanceConstant" -> None,
			"SymbolicRationalHead" -> R,
			"TrialSolutionCoefficient" -> A
		|>;
		
		(* standardise the equations *)
		state["Equation"] = ((# /. Equal -> Subtract) == 0) & /@ ieqns;
		
		(* Provision for custom wave constants *)
		If[
			state["Options"] @ "WaveConstants" =!= Automatic,
				If[Length[state["Options"] @ "WaveConstants"] != Length[state @ "Variables"],
					Message[GERFSolve::ConstantsLengthMismatch];
					Throw @ $Failed];
				state["WaveConstant"] = AssociationThread[
					state @ "Variables"  -> state["Options"] @ "WaveConstants"]];
		
		(* top level *)
		Catch[
			(* convert the eqns to ODE using wave transformation and update state: *)
			state["ODE"] = ReducetoODE[state];
			(* calculate balance constant *)
			state["BalanceConstant"] = AssociationThread[
				Range @ state["Length"] ->
				If[state["Options"] @ "BalanceConstant" === Automatic,
					BalanceConstant[ieqns, us],
					state["Options"] @ "BalanceConstant"]];
			(* and validate it *)
			If[
				!MatchQ[Values @ state["BalanceConstant"], {_Integer ..}] ||
					!AllTrue[Values @ state["BalanceConstant"], Positive],
				Message[GERFSolve::GERFPackageError, "Valid balance constants could not be calculated. Consider providing the values using \"BalanceConstant\" option."];
				Throw[$Failed]];
			(* update state with the trial solutions as functions of wave transform variable, eta *)
			state["TrialSolution"] = AssociationMap[
				TrialSolution[#, state]&, Range @ state["Length"]];
			(* make the auxiliary polynomial which is to be ultimately solved *)
			state["AuxiliaryPolynomial"] = AuxiliaryPolynomial[state];
			(* and finally solve it *)
			SolveAuxiliaryPolynomial[state]]]
GERFSolve[eqn_, u_[vars__], opts : OptionsPattern[]] :=
	GERFSolve[{eqn}, {u[vars]}, opts]


ClearAll[BalanceConstant];
BalanceConstant::usage =
	"BalanceConstant[eqn, u[x, t]] calculates the balancing constant of the given eqn in u[x, t].
BalanceConstant[{eqn1, eqn2, ..., eqnk}, {u1[x, t], ..., uk[x, t]}] calculates the balancing constants of the given equations {eqn1, eqn2, ..., eqnj} in u1[x, t], ..., uk[x, t]}.";
BalanceConstant::GERFPackageError =
	"`1`";


BalanceConstant[ieqns_List, us : {(_[vars__]) ..}] :=
	Module[
		{sys = {}, k, eqns, m, sols, ld, nld, expr, terms, cand, funcs},
		
		funcs = Head /@ us;
		eqns = ((# /. Equal -> Subtract) == 0) & /@ ieqns;
		
		For[k = 1, k <= Length[ieqns], k++,
			expr = Expand @ If[Head[eqns[[k]]] === Equal,
				Subtract @@ eqns[[k]], 
				eqns[[k]]];
			terms = If[Head[expr] === Plus, List @@ expr, {expr}];
			
			(* get degrees  *)
			ld = Simplify[GetDegreeofTerm[#, m, funcs, {vars}] & /@ Select[terms, LinearQ[#, funcs, {vars}]&]];
			nld = Simplify[GetDegreeofTerm[#, m, funcs, {vars}] & /@ Select[terms, Not @* (LinearQ[#, funcs, {vars}]&)]];
			
			If[Length[ld] == 0 || Length[nld] == 0,
				Message[BalanceConstant::GERFPackageError, "Balance constant could not be calculated for equation " <> ToString[k] <> "."];
				Return[$Failed]];

			AppendTo[
				sys, 
				Last[SortBy[ld, # /. {m[_] :> 100, _Symbol :> 1} &]] ==
					Last[SortBy[nld, # /. {m[_] :> 100, _Symbol :> 1} &]]]];
		
		(* solve the simultaneous system for all m[k] *)
		sols = Solve[sys, Table[m[i], {i, 1, Length[ieqns]}]];
		
		If[Length[sols] > 0,
			(* and grab the largest solution *)
			cand = TakeLargestBy[sols, Length, 1][[1]];
			If[Length[cand] == Length[ieqns],
				Last /@ cand,
				Message[BalanceConstant::GERFPackageError, "No valid balance constants found for the system."];
				Return[$Failed]]]]
BalanceConstant[eqn_, u_[vars__]] :=
	Tr @ BalanceConstant[{eqn}, {u[vars]}]


(* ::Section:: *)
(*Package Footer*)


End[];
