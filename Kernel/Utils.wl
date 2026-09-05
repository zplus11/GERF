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


GetDegreeofTerm[term_, m_, funcs_, vars_] :=
	Module[{j, subbed, drules},
		
		(* substitute all functions u[k] with K^m[k] *)
		drules = Flatten @ Table[With[{k = i}, {
			funcs[[k]] @@ vars -> Power[j, m[k]],
			Derivative[orders__][funcs[[k]]] @@ vars :> Power[j, m[k] + Total[{orders}]],
			
			(* a standard fractional derivative of order alpha transforms 
			   to a first-order ODE derivative, so its degree is m[k] + 1 *)
			(FractionalD | CaputoD)[funcs[[k]] @@ vars, {_, alpha_}] :> Power[j, m[k] + 1]
		}], {i, 1, Length[funcs]}];
		
		subbed = term /. drules;
		Simplify[Exponent[subbed, j]]]


FractionalPDEQ[eqn_] := Not @ FreeQ[eqn, FractionalD | CaputoD]


FractionalOrders[eqn_] :=
	Return @ DeleteDuplicates @
		Cases[eqn, (FractionalD | CaputoD)[_, {_, alpha_}] :> alpha, Infinity]


(* ::Text:: *)
(*For extracting balance constants:*)


LinearQ[term_, funcs_, vars_] :=
	Module[{j, subbed, lrules},
		
		(* Substitute all functions u[k] with K to check for linearity (degree 1) *)
		lrules = Flatten @ Table[{
			funcs[[k]] @@ vars -> j,
			Derivative[__][funcs[[k]]] @@ vars :> j,
			(FractionalD | CaputoD)[funcs[[k]] @@ vars, __] :> j
		}, {k, 1, Length[funcs]}];
		
		subbed = term /. lrules;
		Simplify[Exponent[subbed, j]] == 1]


(* ::Text:: *)
(*Integrate the equation:*)


IntegrateEquation[eqn_, state_] :=
	Module[
		{lhs},
		
		lhs = First @ eqn;
		
		While[
			Head @ lhs != Integrate,
			lhs = Integrate[lhs, state @ "Eta"]];
		
		Return[lhs == 0]]


(* ::Text:: *)
(*Reducing the PDE to ODE:*)


ReducetoODE[state_] :=
	Module[
		{dr, interim},
		
		(* applying the derivative rule given as: *)
			dr = Flatten @ Table[
				With[{k = i}, {(FractionalD | CaputoD)[state["Function"][k]@@state["Variables"], {var_, alpha_}] :>
					state["WaveConstant"] @ var * Derivative[1][state["AuxiliaryFunction"][k]][state["Eta"]],
					(* only one FractionalD at a time please *)
					(* D^alpha u_x = x^(1-alpha) U'(eta) deta/dx
						= x^(1-alpha) U'(eta) a x^(alpha-1)
						= a U'(eta) *)
					Derivative[orders__][state["Function"][k]]@@state["Variables"] :>
						Times @@ Power[state["WaveConstant"] /@ state["Variables"], {orders}] *
							Derivative[Total[{orders}]][state["AuxiliaryFunction"][k]][state["Eta"]]}],
					(* simply using this rule, for example u_x converts into 
						a*U', or u_xxy -> a^2bU''', and so on *)
					{i, 1, state @ "Length"}];
		(* to the equation and return *)
			interim = state["Equation"] /. dr /. Table[With[{k = i},
				state["Function"][k]@@state["Variables"] :> state["AuxiliaryFunction"][k][state["Eta"]]],
				{i, 1, state @ "Length"}];
			
			IntegrateEquation[#, state] & /@ interim]


(* ::Text:: *)
(*Constructing the ansatz for Subscript[U, k]:*)


TrialSolution[k_, state_] :=
	Module[
		{A, R, eta},
		
		A = state @ "TrialSolutionCoefficient";
		R = state @ "SymbolicRationalHead";
		
		(* a function of eta *)
		eta = state @ "Eta";
		Return @ Function[
			eta,
			Sum[
				A[i, k] R[eta]^i, 
				{i, -state["BalanceConstant"] @ k, state["BalanceConstant"] @ k}]]]


(* ::Text:: *)
(*Construct the auxiliary polynomials:*)


AuxiliaryPolynomial[state_] :=
	ExpandAll[
		ReplaceAll[state @ "ODE",
			Table[With[{k = i},
			state["AuxiliaryFunction"][k] :> state["TrialSolution"][k]],
			{i, 1, state @ "Length"}]]]


(* ::Text:: *)
(*Algebraic system solver*)


SolveAuxiliaryPolynomial[state_] :=
	Module[
		{tosolve, sys, u, v, sols, vars, forms, pairs},

		tosolve = ReplaceAll[
			state["AuxiliaryPolynomial"],
			state["SymbolicRationalHead"] ->
				state["Options"]["RationalFunction"]];

		sys = Thread[
			CoefficientList[
				Expand @ Numerator @ Together @ TrigToExp @
					(First /@ tosolve) /. {
						Exp[d_. * state["Eta"]] :>
							u^Re[d] v^Im[d],
						state["Eta"] :> u},
				{u, v}] == 0];

		vars = Flatten @ Table[
			state["TrialSolutionCoefficient"][i, k],
			{k, state["Length"]},
			{i,
				-state["BalanceConstant"][k],
				 state["BalanceConstant"][k]}];

		sols = Select[
			Solve @ Reduce[sys, vars],
			Length @ DeleteDuplicates @ #[[All, 2]] > 1 &];

		forms = Table[
			Table[
				state["TrialSolution"][k][FormWaveTransformation[state]],
				{k, state["Length"]}] /.
					state["SymbolicRationalHead"] ->
						state["Options"]["RationalFunction"] /.
							sol,
				{sol, sols}];

		pairs = Select[
			Transpose[{sols, forms}],
				!FreeQ[Last[#], Alternatives @@ state["Variables"]] &&
				FreeQ[Last[#],
					Indeterminate | ComplexInfinity |
					Infinity | DirectedInfinity] &];

		If[pairs === {}, Throw[{}]];

		{sols, forms} = Transpose @ pairs;

		pairs = MapThread[{#1,
				Thread[
					Through[
						(state["Function"] /@
							Range[state["Length"]]) @@
							state["Variables"]
					] -> #2]} &,
			{sols, forms}];

		Throw @ CleanSymbols[
			Switch[
				state["Options"]["OutputMode"],
				"SolutionSets", First /@ pairs,
				All, Transpose[{First /@ pairs, Last /@ pairs}],
				_, Last /@ pairs],
			state]]


(* ::Text:: *)
(*Make solution forms*)


FormWaveTransformation[state_] :=
	Module[
		{raw, fds, defaulters},
		
		(* all cases of fractional orders: *)
		raw = DeleteDuplicates @ Cases[state @ "Equation",
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
		Sum[state["WaveConstant"][k] * k^fds[k] / fds[k] , {k, Keys @ fds}] +
			Sum[state["WaveConstant"][k] * k, {k, Complement[state @ "Variables", Keys @ fds]}]]


CleanSymbols[expr_, state_] :=
	Module[
		{\[ScriptCapitalA], \[ScriptW], A, w},
		A = state["TrialSolutionCoefficient"];
		w = state["WCH"];
		
		expr /. {A[n_, k_] :> 
			If[state["Length"] == 1, Subscript[\[ScriptCapitalA], n], Subsuperscript[\[ScriptCapitalA], n, state["Function"][k]]],
			w[x_] :> Subscript[\[ScriptW], x]} /.
			sym_Symbol /; StringMatchQ[Context[sym], "*Private*"] :> 
				Symbol[StringSplit[SymbolName[sym], "$"][[1]]]]


(* ::Section:: *)
(*Package Footer*)


End[];
