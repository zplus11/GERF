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
		{w, A, R}, (* initialise A,R here so each trial solution has same instance of R *)
		
		$GERFState = <|
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
			"AuxiliaryFunction" -> None,
			"Eta" -> None,
			"WaveConstant" -> AssociationMap[w, {vars}],
			"WCH" -> w, (* ad hoc provision *)
			"BalanceConstant" -> None,
			"SymbolicRationalHead" -> R,
			"TrialSolutionCoefficient" -> A
		|>;
		
		(* standardise the equations *)
		$GERFState @ "Equation" = ((# /. Equal -> Subtract) == 0) & /@ ieqns;
		
		(* Provision for custom wave constants *)
		If[
			$GERFState["Options"] @ "WaveConstants" =!= Automatic,
				If[Length[$GERFState["Options"] @ "WaveConstants"] != Length[$GERFState @ "Variables"],
					Message[GERFSolve::ConstantsLengthMismatch];
					Throw @ $Failed];
				$GERFState @ "WaveConstant" = AssociationThread[
					$GERFState @ "Variables"  -> $GERFState["Options"] @ "WaveConstants"]];
		
		(* top level *)
		Catch[
			(* convert the eqns to ODE using wave transformation and update state: *)
			$GERFState @ "ODE" = ReducetoODE[];
			(* calculate balance constant *)
			$GERFState @ "BalanceConstant" = AssociationThread[
				Range @ $GERFState["Length"] ->
				If[$GERFState["Options"] @ "BalanceConstant" === Automatic,
					BalanceConstant[],
					$GERFState["Options"] @ "BalanceConstant"]];
			(* and validate it *)
			If[
				!MatchQ[Values @ $GERFState["BalanceConstant"], {_Integer ..}] ||
					!AllTrue[Values @ $GERFState["BalanceConstant"], Positive],
				Message[GERFSolve::GERFPackageError, "Valid balance constants could not be calculated. Consider providing the values using \"BalanceConstant\" option."];
				Throw[$Failed]];
			(* update state with the trial solutions as functions of wave transform variable, eta *)
			$GERFState @ "TrialSolution" = AssociationMap[
				TrialSolution, Range @ $GERFState @ "Length"];
			(* make the auxiliary polynomial which is to be ultimately solved *)
			$GERFState @ "AuxiliaryPolynomial" = AuxiliaryPolynomial[];
			(* and finally solve it *)
			SolveAuxiliaryPolynomial[]]]

GERFSolve[eqn_, u_[vars__], opts : OptionsPattern[]] :=
	GERFSolve[{eqn}, {u[vars]}, opts]


(* ::Section:: *)
(*Package Footer*)


End[];
