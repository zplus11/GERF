(* ::Package:: *)

(* ::Section:: *)
(*Package Header*)


BeginPackage["Taggar`GERF`"];


(* ::Text:: *)
(*Declare your public symbols here:*)


GERFSolve::usage =
	"GERFSolve[eqn, u[x, t]] solves the given eqn in u[x, t] using GERF expansion technique.
GERFSolve[{eqn1, eqn2, ..., eqnk}, {u1[x, t], ..., uk[x, t]}] solves the given {eqn1, eqn2, ..., eqnj} in u1[x, t], ..., uk[x, t]} using GERF expansion technique.";


Begin["`Private`"];


(* ::Section:: *)
(*Definitions*)


(* ::Text:: *)
(*Define your public and private symbols here:*)


Get /@ {
	"Taggar`GERF`Utils`",
	"Taggar`GERF`Solver`"
}


GERFSolve::GERFPackageError = "`1`";
GERFSolve::InvalidFractionalDerivatives = "Multiple fractional orders received for these dimensions: `1`. This is not allowed.";
GERFSolve::ConstantsLengthMismatch = "The number of wave constants provided in the \"WaveConstants\" option does not match the number of independent variables in the equation(s).";


(* ::Section::Closed:: *)
(*Package Footer*)


End[];
EndPackage[];
