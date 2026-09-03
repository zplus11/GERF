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
			"WCH" -> w, (* ad hoc provision *)
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
					BalanceConstant[state],
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


(* ::Section:: *)
(*Package Footer*)


End[];
