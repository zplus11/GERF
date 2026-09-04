(* ::Package:: *)

(* ::Section:: *)
(*Package Header*)


Begin["Taggar`GERF`Private`"];


(* ::Section:: *)
(*Definitions*)


(* ::Text:: *)
(*For initial checks of input equation:*)


NonlinearQ[eqn_, u_[vars__]] :=
	Module[
		{expr, pat},
		
		(* has derivative terms or not: *)
		If[Not[Length[{vars}] > 1 &&
			Not @ FreeQ[eqn, Derivative | FractionalD | CaputoD]],
				Throw @ False];
		(* lhs == rhs -> lhs-rhs == 0 *)
		expr = Expand @ If[Head[eqn] === Equal, Subtract @@ eqn, eqn];
		(* testing the pattern *) pat = u[vars] | Derivative[__][u][vars] |
			FractionalD[__] | CaputoD[__];
		(* against expr *)
		Throw @ Not @ FreeQ[expr, Times[a_, b_] /;
			MatchQ[a, pat] && MatchQ[b, pat]]]


GetDegreeofTerm[term_, m_] :=
	Module[{j, subbed, drules},
		
		(* substitute all functions u[k] with K^m[k] *)
		drules = Flatten @ Table[With[{k = i}, {
			$GERFState["Function"][k] @@ $GERFState["Variables"] -> Power[j, m[k]],
			Derivative[orders__][$GERFState["Function"][k]] @@ $GERFState["Variables"] :> Power[j, m[k] + Total[{orders}]],
			
            (* a standard fractional derivative of order alpha transforms 
               to a first-order ODE derivative, so its degree is m[k] + 1 *)
			(FractionalD | CaputoD)[$GERFState["Function"][k] @@ $GERFState["Variables"], {_, alpha_}] :> Power[j, m[k] + 1]
		}], {i, 1, $GERFState["Length"]}];
		
		subbed = term /. drules;
		Simplify[Exponent[subbed, j]]
    ]


FractionalPDEQ[eqn_] := Not @ FreeQ[eqn, FractionalD | CaputoD]


FractionalOrders[eqn_] :=
	Return @ DeleteDuplicates @
		Cases[eqn, (FractionalD | CaputoD)[_, {_, alpha_}] :> alpha, Infinity]


(* ::Text:: *)
(*For extracting balance constants:*)


LinearQ[term_] :=
	Module[{j, subbed, lrules},
		
		(* Substitute all functions u[k] with K to check for linearity (degree 1) *)
		lrules = Flatten @ Table[{
			$GERFState["Function"][k] @@ $GERFState["Variables"] -> j,
			Derivative[__][$GERFState["Function"][k]] @@ $GERFState["Variables"] :> j,
			(FractionalD | CaputoD)[$GERFState["Function"][k] @@ $GERFState["Variables"], __] :> j
		}, {k, 1, $GERFState["Length"]}];
		
		subbed = term /. lrules;
		Simplify[Exponent[subbed, j]] == 1]


(* ::Text:: *)
(*Integrate the equation:*)


IntegrateEquation[eqn_] :=
	Module[
		{lhs},
		
		lhs = First @ eqn;
		
		While[
			Head @ lhs != Integrate,
			lhs = Integrate[lhs, $GERFState @ "Eta"]];
		
		Return[lhs == 0]]


(* ::Text:: *)
(*Function to extract the balance constant of equation when "BalanceConstant" option is set to Automatic:*)


BalanceConstant[] :=
	Module[
		{sys = {}, k, m, sols, ld, nld, expr, terms, cand},
		
		For[k = 1, k <= $GERFState["Length"], k++,
			expr = Expand @ If[Head[$GERFState["Equation"][[k]]] === Equal,
				Subtract @@ $GERFState["Equation"][[k]], 
				$GERFState["Equation"][[k]]];
			terms = If[Head[expr] === Plus, List @@ expr, {expr}];
			
			(* get degrees  *)
			ld = Simplify[GetDegreeofTerm[#, m] & /@ Select[terms, LinearQ]];
			nld = Simplify[GetDegreeofTerm[#, m] & /@ Select[terms, Not @* LinearQ]];
			
			If[Length[ld] == 0 || Length[nld] == 0,
				Message[GERFSolve::GERFPackageError, "Balance constant could not be calculated for equation " <> ToString[k] <> "."];
				Throw[$Failed] (* to top level *)];

			AppendTo[
				sys, 
				Last[SortBy[ld, # /. {m[_] :> 100, _Symbol :> 1} &]] ==
					Last[SortBy[nld, # /. {m[_] :> 100, _Symbol :> 1} &]]]];
		
		(* solve the simultaneous system for all m[k] *)
		sols = Solve[sys, Table[m[i], {i, 1, $GERFState["Length"]}]];
		
		If[Length[sols] > 0,
			(* and grab the largest solution *)
			cand = TakeLargestBy[sols, Length, 1][[1]];
			If[Length[cand] == $GERFState["Length"],
				Last /@ cand,
				Message[GERFSolve::GERFPackageError, "No valid balance constants found for the system."];
				Throw[$Failed] (* to top level *)]]]


(* ::Text:: *)
(*Reducing the PDE to ODE:*)


ReducetoODE[] :=
	Module[
		{dr, U, eta, interim},
		
		$GERFState @ "AuxiliaryFunction" = AssociationMap[U, Range @ $GERFState @ "Length"];
		$GERFState @ "Eta" = eta;
		
		(* applying the derivative rule given as: *)
			dr = Flatten @ Table[
				With[{k = i}, {(FractionalD | CaputoD)[$GERFState["Function"][k]@@$GERFState["Variables"], {var_, alpha_}] :>
					$GERFState["WaveConstant"] @ var * Derivative[1][U[k]][eta],
					(* only one FractionalD at a time please *)
					(* D^alpha u_x = x^(1-alpha) U'(eta) deta/dx
						= x^(1-alpha) U'(eta) a x^(alpha-1)
						= a U'(eta) *)
					Derivative[orders__][$GERFState["Function"][k]]@@$GERFState["Variables"] :>
						Times @@ Power[$GERFState["WaveConstant"] /@ $GERFState["Variables"], {orders}] *
							Derivative[Total[{orders}]][U[k]][eta]}],
					(* simply using this rule, for example u_x converts into 
						a*U', or u_xxy -> a^2bU''', and so on *)
					{i, 1, $GERFState @ "Length"}];
		(* to the equation and return *)
			interim = $GERFState["Equation"] /. dr /. Table[With[{k = i},
				$GERFState["Function"][k]@@$GERFState["Variables"] :> U[k][eta]],
				{i, 1, $GERFState @ "Length"}];
			
			IntegrateEquation /@ interim]


(* ::Text:: *)
(*Constructing the ansatz for Subscript[U, k]:*)


TrialSolution[k_] :=
	Module[
		{A, R, eta},
		
		A = $GERFState @ "TrialSolutionCoefficient";
		R = $GERFState @ "SymbolicRationalHead";
		
		(* a function of eta *)
		eta = $GERFState @ "Eta";
		Return @ Function[
			eta,
			Sum[
				A[i, k] R[eta]^i, 
				{i, -$GERFState["BalanceConstant"] @ k, $GERFState["BalanceConstant"] @ k}]]]


(* ::Text:: *)
(*Construct the auxiliary polynomials:*)


AuxiliaryPolynomial[] :=
	ExpandAll[
		ReplaceAll[$GERFState @ "ODE",
			Table[With[{k = i},
			$GERFState["AuxiliaryFunction"][k] :> $GERFState["TrialSolution"][k]],
			{i, 1, $GERFState @ "Length"}]]]


(* ::Text:: *)
(*Algebraic system solver*)


SolveAuxiliaryPolynomial[] :=
	Module[
		{tosolve, sys, u, v, sols, vars, forms, pairs},

		tosolve = ReplaceAll[
			$GERFState["AuxiliaryPolynomial"],
			$GERFState["SymbolicRationalHead"] ->
				$GERFState["Options"]["RationalFunction"]];

		sys = Thread[
			CoefficientList[
				Expand @ Numerator @ Together @ TrigToExp @
					(First /@ tosolve) /. {
						Exp[d_. * $GERFState["Eta"]] :>
							u^Re[d] v^Im[d],
						$GERFState["Eta"] :> u},
				{u, v}] == 0];

		vars = Flatten @ Table[
			$GERFState["TrialSolutionCoefficient"][i, k],
			{k, $GERFState["Length"]},
			{i,
				-$GERFState["BalanceConstant"][k],
				 $GERFState["BalanceConstant"][k]}];

		sols = Select[
			Solve @ Reduce[sys, vars],
			Length @ DeleteDuplicates @ #[[All, 2]] > 1 &];

		forms = Table[
			Table[
				$GERFState["TrialSolution"][k][FormWaveTransformation[]],
				{k, $GERFState["Length"]}] /.
					$GERFState["SymbolicRationalHead"] ->
						$GERFState["Options"]["RationalFunction"] /.
							sol,
				{sol, sols}];

		pairs = Select[
			Transpose[{sols, forms}],
				!FreeQ[Last[#], Alternatives @@ $GERFState["Variables"]] &&
				FreeQ[Last[#],
					Indeterminate | ComplexInfinity |
					Infinity | DirectedInfinity] &];

		If[pairs === {}, Throw[{}]];

		{sols, forms} = Transpose @ pairs;

		pairs = MapThread[{#1,
				Thread[
					Through[
						($GERFState["Function"] /@
							Range[$GERFState["Length"]]) @@
							$GERFState["Variables"]
					] -> #2]} &,
			{sols, forms}];

		Throw @ CleanSymbols[
			Switch[
				$GERFState["Options"]["OutputMode"],
				"SolutionSets", First /@ pairs,
				All, Transpose[{First /@ pairs, Last /@ pairs}],
				_, Last /@ pairs],
			$GERFState["TrialSolutionCoefficient"],
			$GERFState["WCH"]]]


(* ::Text:: *)
(*Make solution forms*)


FormWaveTransformation[] :=
	Module[
		{raw, fds, defaulters},
		
		(* all cases of fractional orders: *)
		raw = DeleteDuplicates @ Cases[$GERFState @ "Equation",
			(CaputoD | FractionalD)[_, {x_, alpha_}] :> (x -> alpha), 
			Infinity];
		
		(* if at all, fractional orders must appear uniquely for
			each dimension. hence, the following are troublesome: *)
		defaulters = Select[
			GroupBy[raw, First -> Last, DeleteDuplicates],
				Length @ # > 1 &];
		If[Length @ defaulters > 1,
			Message[GERFSolve::InvalidFractionalDerivatives, defaulters];
			Throw @ $Failed]; (* throw to top level *)
		
		(* finally, when all is well, we now formulate the transformation *)
		fds = Association[raw];
		Sum[$GERFState["WaveConstant"][k] * k^fds[k] / fds[k] , {k, Keys @ fds}] +
			Sum[$GERFState["WaveConstant"][k] * k, {k, Complement[$GERFState @ "Variables", Keys @ fds]}]]


CleanSymbols[expr_, A_, w_] :=
	Module[
		{\[ScriptCapitalA], \[ScriptW]},
		expr /. {A[n_, k_] :> 
			If[$GERFState["Length"] == 1, Subscript[\[ScriptCapitalA], n], Subsuperscript[\[ScriptCapitalA], n, k]],
			w[x_] :> Subscript[\[ScriptW], x]} /.
			sym_Symbol /; StringMatchQ[Context[sym], "*Private*"] :> 
				Symbol[StringSplit[SymbolName[sym], "$"][[1]]]]


(* ::Section:: *)
(*Package Footer*)


End[];
