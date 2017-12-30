(* ::Package:: *)

(* ::Title:: *)
(*Quantum Computing Package*)


(* ::Input:: *)
(*(* QuantumComputing`*)
(*   Jacob Austin Marks*)
(*  Wolfram Research*)
(**)*)


BeginPackage["QuantumComputing`"];


(* ::Section:: *)
(*Usage Messages*)


QuantumDiscreteState::usage =
	"QuantumDiscreteState[{\!\(\*SubscriptBox[
StyleBox[\"qob\", \"TI\"], 
StyleBox[\"1\", \"TR\"]]\) \[Rule] \!\(\*SubscriptBox[
StyleBox[\"vals\", \"TI\"], 
StyleBox[\"1\", \"TR\"]]\), \!\(\*SubscriptBox[
StyleBox[\"qob\", \"TI\"], 
StyleBox[\"2\", \"TR\"]]\) \[Rule] \!\(\*SubscriptBox[
StyleBox[\"vals\", \"TI\"], 
StyleBox[\"2\", \"TR\"]]\)}] yields the discrete "<>
	"quantum state in which quantum objects \!\(\*SubscriptBox[
StyleBox[\"qob\", \"TI\"], 
StyleBox[\"i\", \"TI\"]]\) assume values \!\(\*SubscriptBox[
StyleBox[\"vals\", \"TI\"], 
StyleBox[\"i\", \"TI\"]]\).\n"<>
	
	"QuantumDiscreteState[{{\!\(\*SubscriptBox[
StyleBox[\"qob\", \"TI\"], 
StyleBox[\"1\", \"TR\"]]\), \!\(\*SubscriptBox[
StyleBox[\"qob\", \"TI\"], 
StyleBox[\"2\", \"TR\"]]\), \!\(\*
StyleBox[\"\[Ellipsis]\", \"TR\"]\)} \[Rule] \!\(\*
StyleBox[\"vals\", \"TI\"]\)}] yields the discrete "<>
	"quantum state in which potentially entangled quantum objects \!\(\*SubscriptBox[
StyleBox[\"qob\", \"TI\"], 
StyleBox[\"i\", \"TI\"]]\) assume values \!\(\*
StyleBox[\"vals\", \"TI\"]\).\n"<>
	
	"QuantumDiscreteState[{\"Product\" \[Rule] {\!\(\*SubscriptBox[
StyleBox[\"qstate\", \"TI\"], 
StyleBox[\"1\", \"TR\"]]\), \!\(\*SubscriptBox[
StyleBox[\"qstate\", \"TI\"], 
StyleBox[\"2\", \"TR\"]]\), \!\(\*
StyleBox[\"\[Ellipsis]\", \"TR\"]\)}}] yields the quantum "<>
	" state formed by taking the tensor product of discrete quantum states \!\(\*SubscriptBox[
StyleBox[\"qstate\", \"TI\"], \(i\)]\).\n"<>
	
	"QuantumDiscreteState[{\"Mixture\" \[Rule] <|{\!\(\*SubscriptBox[
StyleBox[\"qstate\", \"TI\"], 
StyleBox[\"1\", \"TR\"]]\)\[Rule] \!\(\*SubscriptBox[
StyleBox[\"prob\", \"TI\"], 
StyleBox[\"1\", \"TR\"]]\), \!\(\*SubscriptBox[
StyleBox[\"qstate\", \"TI\"], 
StyleBox[\"2\", \"TR\"]]\)\[Rule] \!\(\*SubscriptBox[
StyleBox[\"prob\", \"TI\"], 
StyleBox[\"2\", \"TR\"]]\), \[Ellipsis]}|>}] yields "<>
	"the mixed quantum state formed from a statistical ensemble in which \!\(\*SubscriptBox[
StyleBox[\"qstate\", \"TI\"], 
StyleBox[\"i\", \"TI\"]]\) "<>
	"has weight \!\(\*SubscriptBox[
StyleBox[\"prob\", \"TI\"], 
StyleBox[\"i\", \"TI\"]]\).\n"<>
	
	"QuantumDiscreteState[{\"Trace\" \[Rule] {\!\(\*
StyleBox[\"qstate\", \"TI\"]\),{\!\(\*SubscriptBox[
StyleBox[\"qob\", \"TI\"], 
StyleBox[\"1\", \"TR\"]]\), \!\(\*SubscriptBox[
StyleBox[\"qob\", \"TI\"], 
StyleBox[\"2\", \"TR\"]]\), \!\(\*
StyleBox[\"\[Ellipsis]\", \"TR\"]\)}}}] yields the quantum "<>
	"discrete state resulting from tracing out subsystems \!\(\*SubscriptBox[
StyleBox[\"qob\", \"TI\"], 
StyleBox[\"i\", \"TI\"]]\) from \!\(\*
StyleBox[\"qstate\", \"TI\"]\).";
	
QuantumDiscreteStateQ::usage = 
	"QuantumDiscreteStateQ[\!\(\*
StyleBox[\"expr\", \"TI\"]\)] gives True if the head of \!\(\*
StyleBox[\"expr\", \"TI\"]\) is "<>
	"QuantumDiscreteState, and False otherwise.\n";

QuantumPureStateQ::usage = 
	"QuantumPureStateQ[\!\(\*
StyleBox[\"expr\", \"TI\"]\)] gives True if \!\(\*
StyleBox[\"expr\", \"TI\"]\) is a pure quantum state "<>
	"and False otherwise.\n";
	
QuantumMixedStateQ::usage = 
	"QuantumMixedStateQ[\!\(\*
StyleBox[\"expr\", \"TI\"]\)] gives True if \!\(\*
StyleBox[\"expr\", \"TI\"]\) is a mixed quantum state "<>
	"and False otherwise.\n";
	
QuantumDiscreteOperation::usage = 
	"QuantumDiscreteOperation[\!\(\*
StyleBox[\"{\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[
StyleBox[\"s\", \"TI\"],\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"pec\", \"TI\"]\)\!\(\*
StyleBox[\",\", \"TI\"]\)\!\(\*
StyleBox[\" \", \"TI\"]\)\!\(\*
StyleBox[\"inputs\", \"TI\"]\)}] yields the unitary discrete quantum "<>
	"operation with specification \!\(\*
StyleBox[\"spec\", \"TI\"]\) on the quantum objects \!\(\*
StyleBox[\"inputs\", \"TI\"]\)\!\(\*
StyleBox[\".\", \"TI\"]\)\!\(\*
StyleBox[\" \", \"TI\"]\)\n"<>
	
	"QuantumDiscreteOperation[\"Conditional\" \[Rule] {\!\(\*SubscriptBox[
StyleBox[\"cond\", \"TI\"], 
StyleBox[\"1\", \"TR\"]]\) \[Rule] \!\(\*SubscriptBox[
StyleBox[\"qop\", \"TI\"], 
StyleBox[\"1\", \"TR\"]]\), \!\(\*SubscriptBox[
StyleBox[\"cond\", \"TI\"], 
StyleBox[\"2\", \"TR\"]]\) \[Rule] \!\(\*SubscriptBox[
StyleBox[\"qop\", \"TI\"], 
StyleBox[\"2\", \"TR\"]]\), \[Ellipsis]}] "<>
	"yields the conditional operation given by \!\(\*SubscriptBox[
StyleBox[\"qop\", \"TI\"], \(i\)]\) where \!\(\*SubscriptBox[
StyleBox[\"cond\", \"TI\"], \(i\)]\) is the first "<>
	"condition satisfied."<>
	
	"QuantumDiscreteOperation[\"Projection\" \[Rule] {\!\(\*
StyleBox[\"spec\", \"TI\"]\)\!\(\*
StyleBox[\",\", \"TI\"]\) \!\(\*
StyleBox[\"qobs\", \"TI\"]\)}] yields the projective "<>
	"measurement specified by \!\(\*
StyleBox[\"spec\", \"TI\"]\) on quantum objects \!\(\*
StyleBox[\"qobs\", \"TI\"]\)\!\(\*
StyleBox[\".\", \"TI\"]\)\n";

QuantumDiscreteOperationQ::usage = 
	"QuantumDiscreteOperationQ[\!\(\*
StyleBox[\"expr\", \"TI\"]\)] gives True if the head of \!\(\*
StyleBox[\"expr\", \"TI\"]\) is "<>
	"QuantumDiscreteOperation, and False otherwise.\n";

QuantumBlochPlot::usage = 
	"QuantumBlochPlot[\!\(\*
StyleBox[\"qstate\", \"TI\"]\)] plots the Bloch vector corresponding to "<>
	"single qubit state \!\(\*
StyleBox[\"qstate\", \"TI\"]\) on a Bloch Sphere.\n"<>
	
	"QuantumBlochPlot[{\!\(\*SubscriptBox[
StyleBox[\"qstate\", \"TI\"], 
StyleBox[\"1\", \"TR\"]]\),\!\(\*SubscriptBox[
StyleBox[\"qstate\", \"TI\"], 
StyleBox[\"2\", \"TR\"]]\),\!\(\*
StyleBox[\"\[Ellipsis]\", \"TR\"]\)}] plots the Bloch vectors for "<>
	"single qubit states \!\(\*SubscriptBox[
StyleBox[\"qstate\", \"TI\"], 
StyleBox[\"i\", \"TI\"]]\) all on a single Bloch Sphere.\n";

QuantumStateDistance::usage = 
	"QuantumStateDistance[\!\(\*SubscriptBox[
StyleBox[\"qstate\", \"TI\"], 
StyleBox[\"1\", \"TR\"]]\),\!\(\*SubscriptBox[
StyleBox[\"qstate\", \"TI\"], 
StyleBox[\"2\", \"TR\"]]\), \"\!\(\*
StyleBox[\"Method\", \"TI\"]\)\" \[Rule] \!\(\*
StyleBox[\"mthd\", \"TI\"]\)\!\(\*
StyleBox[\"]\", \"TI\"]\)\!\(\*
StyleBox[\" \", \"TI\"]\)yields the "<>
	"distance between quantum states \!\(\*SubscriptBox[
StyleBox[\"qstate\", \"TI\"], 
StyleBox[\"1\", \"TR\"]]\) and \!\(\*SubscriptBox[
StyleBox[\"qstate\", \"TI\"], 
StyleBox[\"2\", \"TR\"]]\) as calculated "<>
	"according to method \!\(\*
StyleBox[\"mthd\", \"TI\"]\)\n";
	
QuantumEntangledObjectsQ::usage = 
	"QuantumEntangledObjectsQ[\!\(\*
StyleBox[\"qstate\", \"TI\"]\), {\!\(\*SubscriptBox[
StyleBox[\"qobs\", \"TI\"], 
StyleBox[\"1\", \"TR\"]]\), \!\(\*SubscriptBox[
StyleBox[\"qobs\", \"TI\"], 
StyleBox[\"2\", \"TR\"]]\)}\!\(\*
StyleBox[\"]\", \"TI\"]\)\!\(\*
StyleBox[\" \", \"TI\"]\)yields true if pure "<>
	"quantum state \!\(\*
StyleBox[\"qstate\", \"TI\"]\) is entangled across the bipartition {\!\(\*SubscriptBox[
StyleBox[\"qobs\", \"TI\"], 
StyleBox[\"1\", \"TR\"]]\), \!\(\*SubscriptBox[
StyleBox[\"qobs\", \"TI\"], 
StyleBox[\"2\", \"TR\"]]\)}\n";
	
QuantumCircuit::usage = 
	"QuantumCircuit[\!\(\*
StyleBox[\"{\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[
StyleBox[\"qop\", \"TI\"], 
StyleBox[\"1\", \"TR\"]], \"TI\"]\)\!\(\*
StyleBox[\",\", \"TI\"]\)\!\(\*
StyleBox[\" \", \"TI\"]\)\!\(\*
StyleBox[SubscriptBox[
StyleBox[\"qop\", \"TI\"], 
StyleBox[\"2\", \"TR\"]], \"TI\"]\)\!\(\*
StyleBox[\",\", \"TI\"]\)\!\(\*
StyleBox[\" \", \"TI\"]\)\!\(\*
StyleBox[\"\[Ellipsis]\", \"TI\"]\)}] yields the quantum circuit described by "<>
	" the ordered sequence of operations \!\(\*
StyleBox[SubscriptBox[
StyleBox[\"qop\", \"TI\"], 
StyleBox[\"1\", \"TR\"]], \"TI\"]\)\!\(\*
StyleBox[\",\", \"TI\"]\)\!\(\*
StyleBox[\" \", \"TI\"]\)\!\(\*
StyleBox[SubscriptBox[
StyleBox[\"qop\", \"TI\"], 
StyleBox[\"2\", \"TR\"]], \"TI\"]\)\!\(\*
StyleBox[\",\", \"TI\"]\)\!\(\*
StyleBox[\" \", \"TI\"]\)\!\(\*
StyleBox[\"\[Ellipsis]\", \"TI\"]\) \n";
	
QuantumPartialTr::usage = 
	"QuantumPartialTr[\!\(\*
StyleBox[\"qstate\", \"TI\"]\), \!\(\*
StyleBox[\"qobs\", \"TI\"]\)\!\(\*
StyleBox[\"]\", \"TI\"]\)\!\(\*
StyleBox[\"  \", \"TI\"]\)yields the quantum state resulting from "<>
	"tracing out quantum objects \!\(\*
StyleBox[\"qobs\", \"TI\"]\) from quantum state \!\(\*
StyleBox[\"qstate\", \"TI\"]\)\n";


(* ::Chapter:: *)
(*Function Definitions*)


Begin["`Private`"];


(* ::Subchapter:: *)
(*QuantumDiscreteState*)


Clear[QuantumDiscreteState];

QuantumDiscreteState::rules="Argument `1` is not of the form \!\(\*
StyleBox[\"qobs\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)\[Rule] \!\(\*
StyleBox[\"state\",\nFontSlant->\"Italic\"]\)";
QuantumDiscreteState::spec="Argument `1` is not a valid state specification";
QuantumDiscreteState::dimension="Argument `1` is not a valid state specification"<>
	" for quantum objects of dimension `2`";
QuantumDiscreteState::arity="Argument `1` is not a valid state specification for"<>
	" a state of `2` quantum objects";
QuantumDiscreteState::arrayDimensions="Argument `1` does not have dimensions compatible"<>
	" with state vector or density matrix representation for a quantum state comprised"<>
	" of `3` `2`-level quantum objects";

QuantumDiscreteState::BaseDim="Each base of `1` for BasisState must fall in Range[0,`2`]";
QuantumDiscreteState::BaseArity="Argument `1` must be a list of length `2` to be a BasisState"<>
	" for given list of quantum objects";
	
QuantumDiscreteState::herm="Input density matrix `1` must be Hermitian";
QuantumDiscreteState::tr="Input density matrix `1` must have unit trace";

QuantumDiscreteState::duplicateWires="Quantum object labels must be unique";
QuantumDiscreteState::listDims="All quantum states must have the same "<>
	"total dimension";

QuantumDiscreteState::traceInputs="Argument `1` must be a list of length 2"<>
	" whose first element is the quantum state from which to trace, and "<>
	" whose second element contains the quantum objects to be removed";

QuantumDiscreteState::mixtureAssoc="Argument `1` is not an association with"<>
	" states as keys and probabilities as values";
QuantumDiscreteState::mixtureKeys="Keys of argument `1` are not valid quantum"<>
	" discrete states";
QuantumDiscreteState::mixtureVals="Values of argument `1` must be either symbolic"<>
	" or numeric probabilities";
QuantumDiscreteState::emptyMixture="Argument `1` must be a non-empty association";
QuantumDiscreteState::qobjs="Quantum states in mixture described by argument `1` must"<>
	" correspond to the same quantum objects";
QuantumDiscreteState::traceWires="All quantum objects to be traced out, `1`, "<>
	"must be subsystems of `2`";

QuantumDiscreteState::prodObjs="All quantum states `1` in tensor product state "<>
	"must contain distinct quantum objects";
QuantumDiscreteState::prodDims="All quantum states `1` in tensor product state "<>
	"must have the same qudit dimension";
	
(* Options *)
Options[QuantumDiscreteState] = {"QuditDimension" -> 2};

(* Transform input into a normal form *)
QuantumDiscreteState /: QuantumDiscreteState[specs_List, OptionsPattern[]] := 
	qstateHelper[specs, OptionValue["QuditDimension"]];
	
(* Extract Properties of the State *)
qstate_QuantumDiscreteState?QuantumDiscreteStateQ[prop_String?StringQ] := 
	extractQStateProperty[qstate, prop];

(* Front End Fornatting of State *)
QuantumDiscreteState/: MakeBoxes[qs:QuantumDiscreteState[specs_Association, OptionsPattern[]], fmt_]:=
	stateVisualize[qs, fmt];
	
(* Symbolic Replacement Rule for Quantum State *)
QuantumDiscreteState /: ReplaceAll[qs:QuantumDiscreteState[specs_Association, OptionsPattern[]], subs_]:=
	stateReplaceAll[qs, subs];


Clear[qstateHelper];

qstateHelper[specs_, dim_] := With[{type = specs[[1,1]]},
	Which[
		type === "QuantumObjects", buildStateFromGround[specs[[1,2]], dim],
		type === "Product", buildProductState[specs[[1,2]], dim],
		type === "Mixture", buildMixedState[specs[[1,2]], dim],
		type === "Trace", traceOutFromState[specs[[1,2]]],
		True, buildStateFromGround[specs, dim]
		]];


ClearAll[vecKeys, dmKeys, discreteStateKeysQ];

vecKeys := {"QuditDimension", "QuantumObjects", "StateVector"};
dmKeys := {"QuditDimension", "QuantumObjects", "DensityMatrix"};

discreteStateKeysQ[assoc_] := 
	If[
		(Complement[Keys[assoc], vecKeys] === {}) ||
		(Complement[Keys[assoc], dmKeys] === {}),
		True, 
		False];


Clear[validDiscreteStateQ];

validDiscreteStateQ[qstate_QuantumDiscreteState[assoc_Association, OptionsPattern[]]] := 
	If[discreteStateKeysQ[qstate[[1]]], True, False];
	
validDiscreteStateQ[expr___] := False;


Clear[validDiscreteStateListQ];

validDiscreteStateListQ[states_List] := Which[
	!AllTrue[states, validDiscreteStateQ], 
	False,
	!AllTrue[states[[All, 1, "QuditDimension"]], 
		# === states[[1, 1, "QuditDimension"]]&],
	(Message[QuantumDiscreteState::listDims];False),
	True,
	True
	];
	
validDiscreteStateListQ[expr___] := False;


Clear[QuantumDiscreteStateQ];
QuantumDiscreteStateQ[expr_] := If[Head[expr] === QuantumDiscreteState, True, False];


Clear[vecQ];

vecQ[state_QuantumDiscreteState] := MemberQ[Keys[state[[1]]], "StateVector"];
vecQ[expr___] := $Failed;


Clear[matQ];

matQ[state_QuantumDiscreteState] := !vecQ[state];
matQ[expr___] := $Failed;


Clear[vecToMat];

vecToMat[vec_] := SparseArray[ConjugateTranspose[{vec}]].SparseArray[{vec}];


Clear[matToVecHelper];

matToVecHelper[row_, pos_] := If[!(row[[pos]] ===  0), {row, pos}, False]


Clear[matToVec];

matToVec[mat_] := Module[
	{row, pos},
	{row, pos} = 
		SelectFirst[
			MapIndexed[matToVecHelper[#1, #2[[1]]]&, mat], !(False===#)&];
	row/Sqrt[row[[pos]]]
	]


Clear[l1Normalize];

l1Normalize[vals_List] := Map[#/(Total @ vals)&, vals] 


Clear[wiresFormat];
wiresFormat[ws_] := If[ListQ[ws], ws, {ws}];


stateInputArrayQ[arr_] := 
	If[ListQ[arr] && 
		AllTrue[Flatten[arr], NumericQ[#] || Head[#] === Symbol &], 
			True, False];


Clear[singleQubitStateQ];

singleQubitStateQ[qs_QuantumDiscreteState] := If[
	qs["QuditDimension"] == 2 &&
	Length @ qs["QuantumObjects"] == 1,
	True,
	False];

singleQubitStateQ[expr___] := False;


Clear[sparseQ];
sparseQ[arr_] := If[Head[arr] === SparseArray, True, False];


(* ::Section:: *)
(*Symbolic Replacement for Quantum State*)


Clear[stateReplaceAll];

stateReplaceAll[qs_, subs_] := Module[
	{dim, objs, key, arr, assoc},
	dim = qs["QuditDimension"];
	objs = qs["QuantumObjects"];
	key = If[vecQ[qs], "StateVector", "DensityMatrix"];
	arr = qs[key];
	
	arr = If[sparseQ[arr], 
		SparseArray[ArrayRules[arr]/.subs, Dimensions[arr]],
		arr/.subs];
	assoc = <|{
		key -> arr,
		"QuditDimension" -> dim,
		"QuantumObjects" -> objs
		}|>;
	QuantumDiscreteState[assoc, "QuditDimension" -> dim]
	];


(* ::Section:: *)
(*Properties of Quantum State*)


Clear[extractQStateProperty];

(* Basics *)
extractQStateProperty[qs_, "QuditDimension"] := qs[[1, "QuditDimension"]];
extractQStateProperty[qs_, "QuantumObjects"] := qs[[1, "QuantumObjects"]];
extractQStateProperty[qs_, "StateVector"] := sVec[qs];	
extractQStateProperty[qs_, "DensityMatrix"] := dMat[qs];

(* Basis *)
extractQStateProperty[qs_, "BasisLabels"] := basisLabels[qs];
extractQStateProperty[qs_, "BasisStates"] := basisStates[qs];

(* Numerical Properties *)
extractQStateProperty[qs_, "VonNeumannEntropy"] := vnEntropy[qs];
extractQStateProperty[qs_, "Purity"] := purity[qs];

(* State Queries *)
extractQStateProperty[qs_, "PureStateQ"] := QuantumPureStateQ[qs];
extractQStateProperty[qs_, "MixedStateQ"] := QuantumMixedStateQ[qs];

(* Bloch Sphere *)
extractQStateProperty[qs_, "BlochCartesianCoordinates"] := blochCartesianCoords[qs];
extractQStateProperty[qs_, "BlochSphericalCoordinates"] := blochSphericalCoords[qs];
extractQStateProperty[qs_, "BlochPlot"] := QuantumBlochPlot[qs];

(* Visualizations *)
extractQStateProperty[qs_, "Plot"] := plot[qs];

extractQStateProperty[expr___] := $Failed;


(* ::Subsection:: *)
(*State Queries*)


ClearAll[QuantumPureStateQ, QuantumMixedStateQ];

QuantumPureStateQ[qs_QuantumDiscreteState] := 
	If[qs["Purity"] === 1, True, False];
QuantumPureStateQ[expr___] := $Failed;

QuantumMixedStateQ[qs_QuantumDiscreteState] := 
	If[qs["Purity"] === 1, False, True];
QuantumMixedStateQ[expr___] := $Failed;


(* ::Subsection:: *)
(*State Vector and Density Matrix*)


ClearAll[dMat, sVec];

sVec[qs_] := If[QuantumPureStateQ[qs], qs[[1, "StateVector"]], 
		(Message[QuantumDiscreteState::mixture];$Failed)];
dMat[qs_] := If[MemberQ[Keys[qs[[1]]], "DensityMatrix"], 
		qs[[1, "DensityMatrix"]], 
		vecToMat[qs[[1, "StateVector"]]]];


(* ::Subsection:: *)
(*Basis Labels and Basis Vectors*)


ClearAll[basisLabels, basisStates];

basisLabels[qs_] := qs["QuantumObjects"] /. List -> Ket;
basisStates[qs_] := Module[
	{dim, n, kets},
	dim = qs["QuditDimension"];
	n = Length @ qs["QuantumObjects"];
	kets = Tuples[Range[0, dim-1], n];
	kets = Map[ToString, kets, {2}];
	kets = Map[StringJoin, kets];
	Map[Ket, kets]
	];


(* ::Subsection:: *)
(*Numerical Properties*)


ClearAll[vnEntropy, purity];

vnEntropy[qs_] := Module[
	{dm, eigs, entropy},
	If[vecQ[qs], Return[0]];
	dm = qs[[1,"DensityMatrix"]];
	eigs = Select[Eigenvalues[dm],#>0&];
	entropy = - Total@Map[# Log[2,#]&,eigs];
	entropy
	];

purity[qs_] := If[vecQ[qs], 
	1, With[{dm = qs[[1,"DensityMatrix"]]}, Abs[Tr[dm.dm]]]];


(* ::Subsection:: *)
(*Bloch Sphere*)


(* ::Subsubsection:: *)
(*Bloch Coordinates*)


ClearAll[blochPureSphericalCoords, blochPureCartesianCoords];
ClearAll[blochMixedSphericalCoords, blochMixedCartesianCoords];
ClearAll[blochSphericalCoords, blochCartesianCoords];
ClearAll[u, v, w, r, \[Theta], \[Phi]];

blochPureSphericalCoords[qs_] := 
	With[{
		\[Alpha] = qs["StateVector"][[1]],
		\[Beta] = qs["StateVector"][[2]]
		}, 
		{1, 2 ArcCos[Abs @ \[Alpha]], Arg[\[Beta]] - Arg[\[Alpha]]}
		];
		
blochPureCartesianCoords[qs_] := Module[
	{r, \[Theta], \[Phi], u, v, w},
	{r, \[Theta], \[Phi]} = blochPureSphericalCoords[qs];
	u = r Sin[\[Theta]] Cos[\[Phi]];
	v = r Sin[\[Theta]] Sin[\[Phi]];
	w = r Cos[\[Theta]];
	{u, v, w}
	];
	
blochMixedCartesianCoords[qs_] := Module[
	{dm, u, v, w},
	dm = Normal @ qs["DensityMatrix"];
	u = Re @ (dm[[1, 2]]);
	v = Im @ (dm[[2, 1]]);
	w = dm[[1, 1]] - dm[[2, 2]];
	{u, v, w}
	];
	
blochMixedSphericalCoords[qs_] := Module[
	{r, \[Theta], \[Phi], u, v, w},
	{u, v, w} = blochMixedCartesianCoords[qs];
	If[u == 0 && v == 0 && w == 0,
		Return[{0, 0, 0}]];
	{r, \[Theta], \[Phi]} = Map[Abs, 
				CoordinateTransformData[
					"Cartesian" -> "Spherical", 
					"Mapping", 
					{u, v, w}
					]
				];
	{r, \[Theta], \[Phi]}
	];
	
blochSphericalCoords[qs_] := With[{
	vals = If[vecQ[qs], 
		blochPureSphericalCoords[qs],
		blochMixedSphericalCoords[qs]]},
	Thread[{"r", "\[Theta]", "\[Phi]"} -> vals]];
	
blochCartesianCoords[qs_] := With[{
	vals = If[vecQ[qs], 
		blochPureCartesianCoords[qs],
		blochMixedCartesianCoords[qs]]},
	Thread[{"u", "v", "w"} -> vals]];


(* ::Subsubsection:: *)
(*Bloch Visualization*)


ClearAll[referenceStates, greatCircles, blochQubit];

referenceStates := Graphics3D[{
	Opacity[0.4], Sphere[],
	Black, Thick, Opacity[1.],
	Line[{{0,1,0},{0,-1,0}}],
	Line[{{0,0,1},{0,0,-1}}],
	Line[{{1,0,0},{-1,0,0}}],
	Text["|0\[RightAngleBracket]", {0,0,1.3}],
	Text["|1\[RightAngleBracket]", {0,0,-1.3}],
	Text["|R\[RightAngleBracket]", {0,1.3,0}],
	Text["|L\[RightAngleBracket]", {0,-1.3,0}],
	Text["|+\[RightAngleBracket]", {1.6,-0.25,0}],
	Text["|-\[RightAngleBracket]", {-1.5,0.3,0}]},
	Boxed->False,
	Axes->False,
	PlotRange->{{-1.7, 1.7},{-1.7, 1.7},{-1.7, 1.7}}];
	
greatCircles := ParametricPlot3D[
	{{Cos[t], Sin[t], 0}, 
	 {0, Cos[t], Sin[t]}, 
	 {Cos[t], 0, Sin[t]}}, 
	 {t, 0, 2 \[Pi]},
	PlotStyle->ConstantArray[{Black, Thin}, 3], 
	Axes->False, 
	Boxed->False,
	ImageSize->Medium, 
	PlotRange->{{-1.7, 1.7},{-1.7, 1.7},{-1.7, 1.7}}];
	
Options[blochQubit] = {"Color" -> Red};

blochQubit[qs_QuantumDiscreteState, OptionsPattern[]] := With[{
	coords = blochCartesianCoords[qs][[All, 2]],
	color = OptionValue["Color"]},
	Graphics3D[{
		color, 
		Arrowheads[.03], 
		Arrow[Tube[{{0,0,0},coords},.01],
		{0, -.01}]},
		Boxed->False,
		PlotRange->{{-1.7, 1.7},{-1.7, 1.7},{-1.7, 1.7}}
		]
	]


Clear[QuantumBlochPlot];

QuantumBlochPlot::qubit="Argument `1` is not a single qubit state";
QuantumBlochPlot::input="Inputs to QuantumBlochPlot must be either "<>
	"quantum discrete states or rules of the form qstate \[Rule] color";
QuantumBlochPlot::ruleQubit="Argument `1` must be a single qubit state";
QuantumBlochPlot::ruleColor="Argument `1` must be a color";
	
QuantumBlochPlot[qs_QuantumDiscreteState] := 
	QuantumBlochPlot[{qs}];
QuantumBlochPlot[qs_QuantumDiscreteState -> col_?ColorQ] := 
	QuantumBlochPlot[{qs -> col}];
	
QuantumBlochPlot[qstates_List] := Module[
	{qubits},
	If[!AllTrue[qstates, 
		validBlochInputQ[#]&],
		Return[$Failed]];
	qubits = Map[If[QuantumDiscreteStateQ[#],
					blochQubit[#], 
					blochQubit[#[[1]], "Color" -> #[[2]]]
					]&, qstates];
	Show[Join[
		{greatCircles, referenceStates},
		 qubits]]
	];
QuantumBlochPlot[expr___] := $Failed


Clear[validBlochInputQ];

validBlochInputQ[expr_] := Which[
	QuantumDiscreteStateQ[expr] &&
	!singleQubitStateQ[expr],
	(Message[QuantumBlochPlot::qubit,expr];Return[False]),
	Head[expr] != Rule,
	(Message[QuantumBlochPlot::input,expr];Return[False]),
	!singleQubitStateQ[expr[[1]]],
	(Message[QuantumBlochPlot::ruleQubit,expr[[1]]];Return[False]),
	!ColorQ[expr[[2]]],
	(Message[QuantumBlochPlot::ruleColor,expr[[2]]];Return[False]),
	True,
	True];
validBlochInputQ[expr___] := False;


(* ::Subsection:: *)
(*Plotting Amplitude and Phase*)


ClearAll[colorBar, unitSquare];

colorBar := DensityPlot[y,{x,0,0.1},{y,0,1},
	AspectRatio->10,
	PlotRangePadding->0,
	PlotPoints->{2,150},
	MaxRecursion->0,
	Frame->True,
	FrameTicks->{{None,All},{None,None}},
	ColorFunction->ColorData[{"GrayTones","Reverse"}]];
	
unitSquare[{i_,j_,abs_,arg_}] := 
	{ColorData["BlueGreenYellow"][Rescale[Abs[arg],{0,2\[Pi]}]], 
	Polygon[{{i-1,j,abs},{i,j,abs},{i,j+1,abs},{i-1,j+1,abs}}]};


(* ::Text:: *)
(*Amplitudes and Phases:*)


Clear[plot];

plot[qs_QuantumDiscreteState] := 
	If[vecQ[qs], plotVec[qs],plotMat[qs]]


ClearAll[chartVec, plotVec];

chartVec[inputs_, labels_] := BarChart[
	Style[#[[1]],ColorData[{"Rainbow",{0, 2 Pi}}]@#[[2]]]&
		/@inputs,ImageSize->400, LegendAppearance->"Column", 
		ChartLabels->(labels), 
		PlotLabel->Style["Basis State Amplitudes and Phases","Title",14],
		AxesLabel->{None, "Amplitude"}];
		
plotVec[qs_] := Module[
	{labels, vec, amps, phases, chartInput, chart},
	labels = qs["BasisStates"];
	vec = qs["StateVector"];
	amps = Abs /@ vec;
	phases = Map[If[#>=0, #,  \[Pi] - #]&, Arg /@ vec];
	chartInput = Thread[{amps, phases}];
	chart = chartVec[chartInput, labels];
	Show[Legended[chart, 
		BarLegend[{"Rainbow", {0, 6.28}}, 
			LegendLabel->"Relative Phase"]]
		]
	];


ClearAll[chartMat, plotMat];

chartMat[inputs_, labels_] := BarChart3D[
	Map[Style[#[[1]],ColorData[{"Rainbow",{0, 1}}]@#[[2]]]&,inputs,{2}],
		ChartLayout->"Grid",
		ImageSize->400, 
		LegendAppearance->"Column", 
		ChartLabels-> labels,
		PlotLabel->Style["Basis State Amplitudes and Phases","Title",14], 
		Method->{"Canvas"->None}, AxesLabel->{None,None,"Amplitude"}];

plotMat[qs_] := Module[
	{labels, mat, amps, phases, chartInput, chart},
	labels = qs["BasisStates"];
	labels = {labels, labels /. Ket ->  Bra};
	mat = qs["DensityMatrix"];
	
	phases = Map[Arg, mat, {2}];
	phases = Map[If[#>=0, #,  \[Pi] - #] &, phases, {2}];
	amps = Map[Abs, mat, {2}];
	chartInput = MapThread[Thread[{#1,#2}]&,{amps, phases}];
	chart = chartMat[chartInput, labels];
	Show[Legended[chart, 
		BarLegend[{"Rainbow", {0, 6.28}}, 
			LegendLabel->"Relative Phase"]]
		]
	];


(* ::Section:: *)
(*State Front End Formatting*)


Clear[stateVisualize];

stateVisualize[qs_QuantumDiscreteState, fmt_] := 
	BoxForm`ArrangeSummaryBox[QuantumDiscreteState,
		qs,
		stateArrayVisualize[qs],
		stateBaseVisual[qs], 
		stateExpandedVisual[qs], 
		fmt];


Clear[stateArrayVisualize];

stateArrayVisualize[qs_] := Which[
	vecQ[qs] && Dimensions[qs[[1,"StateVector"]]][[1]] < 8, 
	MatrixForm @ Normal @ qs[[1,"StateVector"]],
	vecQ[qs], 
	qs[[1,"StateVector"]],
	True,
	qs[[1,"DensityMatrix"]]]


Clear[stateBaseVisual];
ClearAll[stateObjsVisualize, stateDimVisualize];

stateBaseVisual[qs_] := 
	Join[stateObjsVisualize[qs], stateDimVisualize[qs]];
	
stateObjsVisualize[qs_] := 
	{BoxForm`MakeSummaryItem[
		{"Quantum Objects: ",qs[[1,"QuantumObjects"]]},
		 StandardForm]};
		 
stateDimVisualize[qs_] := With[{
	dim = qs[[1,"QuditDimension"]]},
	If[dim != 2,
		{BoxForm`MakeSummaryItem[
			{"Dimension: ", dim},
			 StandardForm]},
		{}
		]
	]


Clear[stateExpandedVisual];

stateExpandedVisual[qs_] := Join[
	purityVisualize[qs],
	entropyVisualize[qs],
	blochVisualize[qs]
	];
	
purityVisualize[qs_] := 
	{BoxForm`MakeSummaryItem[
		{"Quantum Purity: ",qs["Purity"]},
		 StandardForm]};

entropyVisualize[qs_] := 
	{BoxForm`MakeSummaryItem[
		{"Von Neumann Entropy: ",qs["VonNeumannEntropy"]},
		 StandardForm]};

blochVisualize[qs_] := 
	If[singleQubitStateQ[qs],
		{BoxForm`MakeSummaryItem[
			{"Bloch Coordinates: ", blochSphericalCoords[qs]},
			 StandardForm]},
		{}
	]


(* ::Section:: *)
(*Building Quantum State*)


(* ::Subsection:: *)
(*States from Input Specifications*)


Clear[buildStateFromGround];

buildStateFromGround[objects_List, dim_Integer] := Module[
	{formattedObjects, wires, specs, subStates},
	
	(* Check that all inputs are of the form qobs \[Rule] state *)
	If[! AllTrue[objects, groundRuleQ], Return[$Failed]];
	
	(* Put inputs in normal form *)
	formattedObjects = 
		Map[wiresFormat[#[[1]]] -> groundStateFormat[#[[2]]] &, objects];
	
	(* Make sure all state labels are unique *)
	If[!uniqueStateWiresQ[formattedObjects], Return[$Failed]];
	
	(* Check that all state specifications are valid discrete states *)
	If[! AllTrue[specs, discreteStateQ], Return[$Failed]];
	
	(* Check that dimensions agree with state specifications *)
	If[! AllTrue[formattedObjects[[All, 2]], stateDimQ[#, dim]&], Return[$Failed]];
	
	(* Check that number of inputs agree with state specifications *)
	If[! AllTrue[formattedObjects, stateArityQ[#[[2]], Length @ #[[1]], dim]&], Return[$Failed]];
	
	(* Check that if input is a density matrix, then Hermitian and Unit Trace*)
	If[! AllTrue[formattedObjects, notDMQ[#, dim] || (hermitianQ[#[[2]]] && traceIsOneQ[#[[2]]]) &], Return[$Failed]];
		
	subStates = Map[generateSubstate[#[[2]], #[[1]], dim]&, formattedObjects];
	buildProductState[subStates, dim]
	];


Clear[generateSubstate];

generateSubstate[spec_, wires_, dim_] := Module[
	{assoc},
	assoc = <|{
		"QuditDimension" -> dim, 
		"QuantumObjects" -> Sort[wires]}|>;
	If[pureSubstateQ[spec, Length @ wires, dim],
		assoc["StateVector"] = generateSubVec[spec, wires, dim],
		assoc["DensityMatrix"] = generateSubMat[spec, wires, dim]
		];
	QuantumDiscreteState[assoc]
	];


Clear[generateSubVec];

generateSubVec[spec_, wires_, d_] := Module[
	{n, vec, shape, tpLevels},
	n = Length @ wires;
	vec = subVecHelper[spec, n, d];
	shape = ConstantArray[d, n];
	vec = ArrayReshape[vec, shape];
	tpLevels = InversePermutation[Ordering @ wires];
	vec = Normalize @ 
			SparseArray[Flatten[Transpose[vec, tpLevels]]];
	vec
	];


Clear[subVecHelper];

(* Computational Basis States *)
subVecHelper[{"BasisState", bases_}, n_Integer, d_Integer] := basisState[bases, d];

(* Registers *)
subVecHelper["Register", n_Integer, d_Integer] := register[0, n, d];
subVecHelper[{"Register", num_}, n_Integer, d_Integer] := register[num, n, d];


(* States from Input Arrays *)
subVecHelper[arr_List, n_Integer, d_Integer] := arr;

(* +/- Basis States *)
subVecHelper["Plus", 1, 2] := plus;
subVecHelper["Minus", 1, 2] := minus;

(* R/L Basis States *)
subVecHelper["Right", 1, 2] := right;
subVecHelper["Left", 1, 2] := left;

(* Bell States *)
subVecHelper["PsiPlus", 2, 2] := psiPlus;
subVecHelper["PsiMinus", 2, 2] := psiMinus;
subVecHelper["PhiPlus", 2, 2] := phiPlus;
subVecHelper["PhiMinus", 2, 2] := phiMinus;

(* Entangled States *)
subVecHelper["GHZ", n_Integer, 2] := ghzState[n];
subVecHelper["W", n_Integer, 2] := wState[n];

(* Random Pure States *)
subVecHelper["RandomPure", n_Integer, d_Integer] := randPure[n, d];

(* Pure Uniform Superposition *)

subVecHelper["PureUniformSuperposition", n_Integer, d_Integer] :=
	pureSuperpos[n, d];


Clear[generateSubMat];

generateSubMat[spec_, wires_, d_] := Module[
	{n, mat, shape, ord, tpLevels},
	n = Length @ wires;
	mat = subMatHelper[spec, n, d];
	shape = ConstantArray[d, 2 * n];
	mat = ArrayReshape[mat, shape];
	ord = Ordering @ wires;
	tpLevels = InversePermutation[Join[ord, ord + n]];
	mat = Flatten[Transpose[mat, tpLevels]];
	mat = SparseArray @ ArrayReshape[mat, {d^n, d^n}];
	mat
	];


Clear[subMatHelper];

subMatHelper[arr_List, n_Integer, d_Integer] := arr;


Clear[groundStateFormat];

groundStateFormat[{"BasisState", n_Integer}] := {"BasisState", {n}};

groundStateFormat["+"|"plus"|"PLUS"] := "Plus";
groundStateFormat["-"|"minus"|"MINUS"] := "Minus";
groundStateFormat["l"|"L"|"left"|"LEFT"] := "Left";
groundStateFormat["r"|"R"|"right"|"RIGHT"] := "Right";

groundStateFormat["\[Psi]+"|"\[CapitalPsi]+"|"\!\(\*SubscriptBox[\(\[Psi]\), \(+\)]\)"|"\!\(\*SubscriptBox[\(\[CapitalPsi]\), \(+\)]\)"|"psiplus"|"PSIPLUS"] := "PsiPlus";
groundStateFormat["\[Phi]+"|"\[CapitalPhi]+"|"\!\(\*SubscriptBox[\(\[Phi]\), \(+\)]\)"|"\!\(\*SubscriptBox[\(\[CapitalPhi]\), \(+\)]\)"|"phiplus"|"PHIPLUS"] := "PhiPlus";
groundStateFormat["\[Psi]-"|"\[CapitalPsi]-"|"\!\(\*SubscriptBox[\(\[Psi]\), \(-\)]\)"|"\!\(\*SubscriptBox[\(\[CapitalPsi]\), \(-\)]\)"|"psiminus"|"PSIMINUS"] := "PsiMinus";
groundStateFormat["\[Phi]-"|"\[CapitalPhi]-"|"\!\(\*SubscriptBox[\(\[Phi]\), \(-\)]\)"|"\!\(\*SubscriptBox[\(\[CapitalPhi]\), \(-\)]\)"|"phiminus"|"PHIMINUS"] := "PhiMinus";

groundStateFormat["Ghz"|"ghz"] := "GHZ";
groundStateFormat["w"] := "W";

groundStateFormat[expr___] := expr;


Clear[basisState, plus, minus, right, left];
Clear[phiMinus, phiPlus, psiMinus, psiPlus];
Clear[ghzState, wState];
Clear[randPure, pureSuperpos];

plus := {1,1};
minus := {1,-1};

right := {1,I};
left := {1,-I};

psiMinus := {1,0,0,-1};
psiPlus := {1,0,0,1};
phiMinus := {0, 1,-1,0};
phiPlus := {0, 1,1,0};

basisState[bases_, d_] := With[{
	pos = FromDigits[bases, d] + 1,
	n = Length @ bases},
	SparseArray[{pos} -> 1, {d^n}]];
	
ghzState[n_] := 
	SparseArray[{{1} -> 1, {2^n} -> 1},{2^n}];
	
wState[n_] := 
	SparseArray[{i_}/;IntegerQ[Log[2,i-1]]-> 1,{Power[2, n]}];
	
randPure[n_, d_] := RandomComplex[{-1 - I, 1 + I}, d^n];

pureSuperpos[n_, d_] := ConstantArray[1, d^n];

register[num_, n_, d_] := 
	SparseArray[{{num+1} -> 1},{d^n}];


(* ::Subsubsection:: *)
(*Input Spec Queries*)


Clear[groundRuleQ];
groundRuleQ[expr_Rule] := True;
groundRuleQ[expr___] := (Message[QuantumDiscreteState::rules,expr];False);


Clear[uniqueStateWiresQ];

uniqueStateWiresQ[objs_] := If[DuplicateFreeQ[Flatten @ objs[[All, 1]]], 
		True, 
		(Message[QuantumDiscreteState::duplicateWires];False)];


Clear[pureSubstateQ];

pureSubstateQ[arr_?stateInputArrayQ, n_Integer, d_Integer] := With[
	{size = d^n},
	If[Dimensions[arr] === {size, size}, False, True]];
	
pureSubstateQ[expr___] := True;


Clear[discreteStateQ];

discreteStateQ[{"BasisState", bases_}] := 
	If[IntegerQ[bases]||(ListQ[bases]&& AllTrue[bases,IntegerQ]), True, False];

discreteStateQ["Register"] := True;
discreteStateQ[{"Register", num_Integer}] := True;

discreteStateQ["Plus"] := True;
discreteStateQ["Minus"] := True;
discreteStateQ["Left"] := True;
discreteStateQ["Right"] := True;

discreteStateQ["PhiPlus"] := True;
discreteStateQ["PhiMinus"] := True;
discreteStateQ["PsiPlus"] := True;
discreteStateQ["PsiMinus"] := True;

discreteStateQ["GHZ"] := True;
discreteStateQ["W"] := True;

discreteStateQ["RandomPure"] := True;
discreteStateQ["RandomMixed"] := True;

discreteStateQ["PureUniformSuperposition"] := True;
discreteStateQ[arr_?stateInputArrayQ] := True

discreteStateQ[expr___] := (Message[QuantumDiscreteState::spec,expr];False);


Clear[stateDimQ];

stateDimQ[{"BasisState", bases_}, dim_Integer] := 
	If[Max[bases]+1 <= dim, 
		True, 
		(Message[QuantumDiscreteState::BaseDim,bases,(dim-1)];False)];
		
stateDimQ["Register", dim_Integer] := True;
stateDimQ[{"Register", num_Integer}, dim_Integer] := True;

stateDimQ["Plus", 2] := True;
stateDimQ["Minus", 2] := True;
stateDimQ["Left", 2] := True;
stateDimQ["Right", 2] := True;

stateDimQ["PhiPlus", 2] := True;
stateDimQ["PhiMinus", 2] := True;
stateDimQ["PsiPlus", 2] := True;
stateDimQ["PsiMinus", 2] := True;

stateDimQ["GHZ", 2] := True;
stateDimQ["W", 2] := True;

stateDimQ["RandomPure", dim_Integer] := True;
stateDimQ["RandomMixed", dim_Integer] := True;

stateDimQ["PureUniformSuperposition", dim_Integer] := True;
stateDimQ[arr_?stateInputArrayQ, dim_Integer] := True

stateDimQ[expr_, dim_] := (Message[QuantumDiscreteState::dimension,expr,dim];False);


Clear[notDMQ];

notDMQ[spec_, dim_] := With[{
	state = spec[[2]],
	nObjs = Length @ spec[[1]]},
	If[
		!ListQ[state] || 
		!stateInputArrayQ[state] ||
		Dimensions @ state !=  {dim ^ nObjs, dim ^ nObjs},
		True,
		False
		]
		
	];


Clear[hermitianQ];

hermitianQ[arr_] := If[
		HermitianMatrixQ[arr],
		True,
		(Message[QuantumDiscreteState::herm,arr];False)
		]


Clear[traceIsOneQ];

traceIsOneQ[arr_] := 
	If[Tr[arr] === 1,
		True,
		(Message[QuantumDiscreteState::tr,arr];False)];
		
traceIsOneQ[expr___] := True;


Clear[stateArityQ];

stateArityQ[{"BasisState", bases_}, n_, dim_Integer] := 
	If[Length @ bases === n, 
		True, 
		(Message[QuantumDiscreteState::BaseArity,bases,n];False)];

stateArityQ["Register", n_Integer, dim_Integer] := True;
stateArityQ[{"Register", num_Integer}, n_Integer, dim_Integer] := True;

stateArityQ["Plus", 1, 2] := True;
stateArityQ["Minus", 1, 2] := True;
stateArityQ["Left", 1, 2] := True;
stateArityQ["Right", 1, 2] := True;

stateArityQ["PhiPlus", 2, 2] := True;
stateArityQ["PhiMinus", 2, 2] := True;
stateArityQ["PsiPlus", 2, 2] := True;
stateArityQ["PsiMinus", 2, 2] := True;

stateArityQ["GHZ", n_Integer, 2] := True;
stateArityQ["W", n_Integer, 2] := True;

stateArityQ["RandomPure", n_Integer, dim_Integer] := True;
stateArityQ["RandomMixed", n_Integer, dim_Integer] := True;

stateArityQ["PureUniformSuperposition", n_Integer, dim_Integer] := True;

stateArityQ[arr_?stateInputArrayQ, n_Integer, d_Integer] := With[
	{
	 size = d^n,
	 dims = Dimensions[arr]
	 },
	If[dims === ConstantArray[d, n] || dims === {size, size}, 
		True,
		(Message[QuantumDiscreteState::arrayDimensions,arr,d,n];False)]
	];


stateArityQ[expr_, n_, dim_] := (Message[QuantumDiscreteState::arity,expr,n];False);


(* ::Subsection:: *)
(*Mixed States*)


Clear[buildMixedState];

buildMixedState[specs_, dim_Integer] := Module[
	{pairs},
	
	(* Check that input is an association with *)
	(* state keys and prob values *)
	If[!mixtureInputQ[specs], Return[$Failed]];
	
	(* Make sure association isn't empty *)
	If[emptyMixtureQ[specs], Return[$Failed]];
	
	(* If only single component in mixture, return that component *)
	If[Length @ specs === 1, Return[Keys[specs][[1]]]];
	
	(* Make sure all input states describe same quantum objects *)
	If[!sameObjectsQ[specs], Return[$Failed]];
	
	mixture[specs, dim]
	];


Clear[mixture];

mixture[specs_, dim_] := Module[
	{states, probs, ws, reps, dm, assoc},
	states = Keys[specs];
	probs = l1Normalize @ Values[specs];
	ws = states[[1,1]]["QuantumObjects"];
	reps = Map[If[vecQ[#], 
		vecToMat @ #[[1,"StateVector"]],
		#[[1, "DensityMatrix"]]]&, states];
	dm = Total @ MapThread[Times,{reps, probs}];
	assoc = <|{
		"QuditDimension" -> dim,
		"QuantumObjects" -> ws
		}|>;
	If[Tr[dm.dm] === 1, 
		assoc["StateVector"] = matToVec[dm],
		assoc["DensityMatrix"] = dm];
	QuantumDiscreteState[assoc]
	];


(* ::Subsubsection:: *)
(*Mixture Queries*)


Clear[mixtureInputQ];

mixtureInputQ[specs_] := Which[
	!AssociationQ[specs],
	(Message[QuantumDiscreteState::mixtureAssoc,specs];False),
	(*!validDiscreteStateListQ[Keys[specs]],
	(Message[QuantumDiscreteState::mixtureKeys,specs];False),*)
	!AllTrue[Values[specs], Head[#] === Symbol || NumericQ[#] &],
	(Message[QuantumDiscreteState::mixtureVals,specs];False),
	True,
	True];


Clear[emptyMixtureQ];

emptyMixtureQ[assoc_] := 
	If[Length @ assoc === 0,
		(Message[QuantumDiscreteState::emptyMixture,assoc];True),
		False];


Clear[sameObjectsQ];

sameObjectsQ[assoc_Association] := With[{
	qobs = (Keys[assoc])[[All, 1, "QuantumObjects"]]},
	If[AllTrue[qobs, # === qobs[[1]] &],
		True,
		(Message[QuantumDiscreteState::qobjs,assoc];False)]
	]


(* ::Subsection:: *)
(*Product States*)


ClearAll[productStatesObjectOverlap, productStateDims];

productStatesObjectOverlap[states_] := With[
	{objs = Flatten @ Map[#["QuantumObjects"]&, states]},
	If[DuplicateFreeQ[objs], True,
		(Message[QuantumDiscreteState::prodObjs,states];False)
		]
	];
	
productStateDims[states_] := With[
	{dims = Map[#["QuditDimension"]&, states]},
	If[dims === ConstantArray[dims[[1]], Length@states], True,
		(Message[QuantumDiscreteState::prodDims,states];False)
		]
	];


Clear[buildProductState];

buildProductState[subStates_, dim_] := Which[
	!productStatesObjectOverlap[subStates],
	$Failed,
	!productStateDims[subStates],
	$Failed,
	Length @ subStates === 1, subStates[[1]],
	AllTrue[subStates, vecQ], tensorPure[subStates, dim],
	True, tensorMixed[subStates, dim]];


Clear[tensorPure];

tensorPure[subStates_, dim_] := Module[
	{rep, wires, shape, tpLevels, assoc},
	rep = KroneckerProduct @@ subStates[[All, 1, "StateVector"]];
	wires = Flatten @ subStates[[All, 1, "QuantumObjects"]];
	shape = ConstantArray[dim, Length @ wires];
	tpLevels = InversePermutation[Ordering @ wires];
	rep = Flatten @ Transpose[rep, tpLevels];
	assoc = <|{
		"StateVector" -> SparseArray[rep],
		"QuditDimension" -> dim,
		"QuantumObjects" -> Sort[wires]
		}|>;
	QuantumDiscreteState[assoc]
	];


Clear[tensorMixed];

tensorMixed[subStates_, dim_] := Module[
	{nObjs, mats, rep, wires, shape, order, tpLevels, assoc},
	mats = Map[If[vecQ[#], 
		vecToMat @ #[[1,"StateVector"]],
		#[[1, "DensityMatrix"]]]&, subStates];
	rep = KroneckerProduct @@ mats;
	wires = Flatten @ subStates[[All, 1, "QuantumObjects"]];
	nObjs = Length @ wires;
	shape = ConstantArray[dim, 2 * nObjs];
	order = Ordering @ wires;
	tpLevels = InversePermutation[Join[order, order + nObjs]];
	rep = Transpose[rep, tpLevels];
	rep = ArrayReshape[Flatten @ rep, {dim^nObjs, dim^nObjs}];
	assoc = <|{
		"DensityMatrix" -> SparseArray[rep],
		"QuditDimension" -> dim,
		"QuantumObjects" -> Sort[wires]
		}|>;
	QuantumDiscreteState[assoc]
	];


(* ::Subsection:: *)
(*Trace Out Quantum Objects from State*)


Clear[traceOutFromState];

traceOutFromState[specs_] := 
	If[!ListQ[specs] || Length[specs] != 2,
		(Message[QuantumDiscreteState::traceInputs,specs];Return[$Failed]),
		QuantumPartialTr[specs[[1]], specs[[2]]]]


Clear[QuantumPartialTr];

QuantumPartialTr[state_, wires_] := With[{
	ws = If[ListQ[wires], wires, {wires}]},
	Which[
		(*!validDiscreteStateQ[state]*)
		!QuantumDiscreteStateQ[state],
		$Failed,
		!AllTrue[ws, 
			MemberQ[state[[1, "QuantumObjects"]], #]&],
			(Message[QuantumDiscreteState::traceWires,ws,state];$Failed),
		Length @ state["QuantumObjects"] === Length @ wires,
		1,
		True,
		traceHelper[state, ws]
		]
	];

QuantumPartialTr[expr___] := $Failed;


Clear[traceHelper];

traceHelper[state_, {}] := state;

traceHelper[state_, ws_] := Module[
	{newState, dm, dim, wContract, stateWs, pos, n, tpLevels, assoc},
	dm = state["DensityMatrix"];
	dim = state["QuditDimension"];
	wContract = ws[[1]];
	stateWs = state["QuantumObjects"];
	n = Length @ stateWs;
	
	pos = Position[stateWs, wContract][[1,1]];
	dm = ArrayReshape[dm, ConstantArray[dim, 2 * n]];
	tpLevels = InversePermutation[Join[Range[n], Range[n] + n]];
	dm = Transpose[dm, tpLevels];
	dm = TensorContract[dm, {{pos, pos + n}}];
	tpLevels = InversePermutation@Join[Range[n - 1], Range[n - 1] + (n - 1)];
	dm = Transpose[dm, tpLevels];
	dm = ArrayReshape[Flatten @ dm, {dim^(n-1), dim^(n-1)}];
	assoc = <|{
		"QuantumObjects" -> Complement[stateWs, {wContract}],
		"QuditDimension" -> dim
		}|>;
	If[Tr[dm.dm] == 1,
		assoc["StateVector"] = matToVec[dm],
		assoc["DensityMatrix"] = dm];
	newState = QuantumDiscreteState[assoc];
	traceHelper[newState, Drop[ws, 1]]
	]


(* ::Subsection:: *)
(*QuantumStateDistance *)


Clear[QuantumStateDistance];

QuantumStateDistance::inputs="Inputs to QuantumStateDistance must be two "<>
	"quantum discrete states";
QuantumStateDistance::distMeasure="Method `1` is not a valid option for "<>
	"calculating the distance between two states";
QuantumStateDistance::qubit="Argument `1` must be a single qubit state for "<>
	" Euclidean distance method";
QuantumStateDistance::dims="Arguments `1` at position 1 and `2` at position 2 "<>
	"must have the same qudit dimension";
QuantumStateDistance::numObjs="Arguments `1` at position 1 and `2` at position 2 "<>
	"must describe the same number of quantum objects";

(* Options *)
Options[QuantumStateDistance] = {"Method" -> "Fidelity"};

QuantumStateDistance[qs1_QuantumDiscreteState, qs2_QuantumDiscreteState, OptionsPattern[]] := If[
	!sameSizeStatesQ[qs1, qs2], $Failed,
	qdistHelper[qs1, qs2, OptionValue["Method"]]];
	
QuantumStateDistance[expr___] := 
	(Message[QuantumStateDistance::inputs];$Failed);


Clear[sameSizeStatesQ];

sameSizeStatesQ[qs1_, qs2_] := Which[
	qs1["QuditDimension"] != qs2["QuditDimension"],
	(Message[QuantumStateDistance::dims,qs1,qs2];False),
	Length @ qs1["QuantumObjects"] != Length @ qs2["QuantumObjects"],
	(Message[QuantumStateDistance::numObjs,qs1,qs2];False),
	True,
	True]


Clear[qdistHelper];

qdistHelper[qs1_, qs2_, "Fidelity"] := fidelityDist[qs1, qs2];
qdistHelper[qs1_, qs2_, "Trace"] := traceDist[qs1, qs2];
qdistHelper[qs1_, qs2_, "BuresAngle"] := buresAngleDist[qs1, qs2];
qdistHelper[qs1_, qs2_, "HilbertSchmidt"] := hilbertSchmidtDist[qs1, qs2];
qdistHelper[qs1_, qs2_, "Euclidean"] := blochDist[qs1, qs2];
qdistHelper[qs1_, qs2_, method_] := 
	(Message[QuantumStateDistance::distMeasure];$Failed);


Clear[fidelityDist, traceDist, buresAngleDist, hilbertSchmidtDist, blochDist];

fidelityDist[qs1_, qs2_] := With[{
	dm1 = qs1["DensityMatrix"],
	rootDM2 = MatrixPower[qs2["DensityMatrix"], 1/2]},
	Re @ Chop[1 - Tr[MatrixPower[rootDM2. dm1. rootDM2,1/2]], 10^(-6)]
	];

traceDist[qs1_, qs2_] := Module[
	{\[Rho], \[Sigma], eigs},
	\[Rho] = qs1["DensityMatrix"];
	\[Sigma] = qs2["DensityMatrix"];
	eigs = Eigenvalues[\[Rho] - \[Sigma]];
	Re @ Total @ (Abs /@ eigs)/2
	];
	
(* Bures Angle - Statistical Distance *)
buresAngleDist[qs1_, qs2_] := 
	Re @ ArcCos[Sqrt[fidelityDist[qs1, qs2]]];
	
hilbertSchmidtDist[qs1_, qs2_] := Module[
	{\[Rho], \[Sigma]},
	\[Rho] = qs1["DensityMatrix"];
	\[Sigma] = qs2["DensityMatrix"];
	Re @ Sqrt @ Tr @ MatrixPower[\[Rho] - \[Sigma], 2]
	];

blochDist[qs1_, qs2_] := Module[
	{vec1, vec2},
	If[!singleQubitStateQ[qs1], 
		(Message[QuantumStateDistance::qubit,qs1];Return[$Failed])];
	If[!singleQubitStateQ[qs2], 
		(Message[QuantumStateDistance::qubit,qs2];Return[$Failed])];
	vec1 = (qs1["BlochCartesianCoordinates"])[[All, 2]];
	vec2 = (qs2["BlochCartesianCoordinates"])[[All, 2]];
	Re @ EuclideanDistance[vec1, vec2]
	];


(* ::Subsection:: *)
(*QuantumEntangledObjectsQ*)


Clear[QuantumEntangledObjectsQ];

QuantumEntangledObjectsQ::inputs="Arguments to QuantumEntangledObjectsQ must be "<>
	"a pure quantum discrete state qstate and a bipartition of the quantum objects "<>
	"in qstate";
QuantumEntangledObjectsQ::qstate="Argument `1` at position 1 is not a quantum discrete state";
QuantumEntangledObjectsQ::pureState="Argument `1` at position 1 must be pure";
QuantumEntangledObjectsQ::splitList="Argument `1` at position 2 must be a list";
QuantumEntangledObjectsQ::splitLength="Argument `1` at position 2 must have length 2";
QuantumEntangledObjectsQ::qObjs="All elements of argument `2` at position `3` in `4` "<>
	"must be quantum objects in `1`";
QuantumEntangledObjectsQ::unique="All elements of Flatten @ `1` at position 2 must be unique";
QuantumEntangledObjectsQ::bipartition="Argument `2` at position 2 must form a "<>
	"valid bipartition of the quantum objects in `1`";
	
QuantumEntangledObjectsQ[qs_, split_] := Module[
	{redState, tracedObjects, redDM, concurrence},
	
	If[!validEntangledStateQ[qs],Return[$Failed]];
	If[!validBipartitionQ[qs, split],Return[$Failed]];

	tracedObjects = split[[1]];
	redState = QuantumPartialTr[qs, tracedObjects];
	redDM = redState["DensityMatrix"];
	concurrence = Sqrt[2(1 - Tr[redDM.redDM])];
	If[concurrence > 0, True, False]
];


QuantumEntangledObjectsQ[expr___] := 
	(Message[QuantumEntangledObjectsQ::inputs];$Failed);


(* ::Subsubsection:: *)
(*Entangled Object Queries*)


Clear[validEntangledStateQ];

validEntangledStateQ[qs_QuantumDiscreteState] := If[
	qs["PureStateQ"],
	True,
	(Message[QuantumEntangledObjectsQ::pureState,qs];False)];

validEntangledStateQ[expr___] := 
	(Message[QuantumEntangledObjectsQ::qstate,expr];False);


Clear[validBipartitionQ];

validBipartitionQ[qs_, split_List] := Which[
	Length[split] != 2,
	(Message[QuantumEntangledObjectsQ::splitLength,split];False),
	!AllTrue[split[[1]],MemberQ[qs["QuantumObjects"], #]&],
	(Message[QuantumEntangledObjectsQ::qObjs,qs,split[[1]],1,split];False),
	!AllTrue[split[[2]],MemberQ[qs["QuantumObjects"], #]&],
	(Message[QuantumEntangledObjectsQ::qObjs,qs,split[[2]],2,split];False),
	!DuplicateFreeQ[Flatten @ split],
	(Message[QuantumEntangledObjectsQ::unique,split];False),
	!Complement[qs["QuantumObjects"], Flatten@split] === {},
	(Message[QuantumEntangledObjectsQ::bipartition,qs,split];False),
	True,
	True
	];
	
validBipartitionQ[qs_, split_] := 
	(Message[QuantumEntangledObjectsQ::splitList,split];False);


(* ::Chapter:: *)
(*QuantumDiscreteOperation*)


Clear[QuantumDiscreteOperation];

QuantumDiscreteOperation::op="Argument `1` at position 1 is not a valid specification";
QuantumDiscreteOperation::arity="Number of wires not equal to number of gate inputs for `1`";
QuantumDiscreteOperation::wires="Dimensions of argument `1` at position 1 are not compatible with the specified operations.";
QuantumDiscreteOperation::power="Not an accepted power for applying operator at position 1.";
QuantumDiscreteOperation::opType="Argument `1` is not a valid type of quantum operation";
QuantumDiscreteOperation::unitary="Argument `1` is not a unitary matrix";
QuantumDiscreteOperation::unitarySpec="Argument `1` is not a valid specification for a unitary operation";
QuantumDiscreteOperation::listDims="All elements of list must act on qudits of the same dimension";

QuantumDiscreteOperation::opWiresNotInState="Input quantum objects `1` to operation must be "<>
	"elements of state `2`";
QuantumDiscreteOperation::opStateDims="Dimension of state qudits `2 and of operation qudits `1` must be equal";

(* Options *)
Options[QuantumDiscreteOperation] = {"QuditDimension" -> 2, "Power" -> 1};

(* Extract Properties of the Operation *)
qop_QuantumDiscreteOperation?QuantumDiscreteOperationQ[prop_String?StringQ] := extractQopProperty[qop, prop];

(* Transform input into a normal form *)
QuantumDiscreteOperation /: QuantumDiscreteOperation[type_ -> specs_List, OptionsPattern[]] := 
	qopHelper[type, specs, OptionValue["QuditDimension"], OptionValue["Power"]];
	
QuantumDiscreteOperation /: QuantumDiscreteOperation[specs_List, OptionsPattern[]] := 
	qopHelper["Unitary", specs, OptionValue["QuditDimension"], OptionValue["Power"]];
	
(* Exponentiation of Unitary Operations *)	
QuantumDiscreteOperation /: ((qop:QuantumDiscreteOperation[specs_Association, OptionsPattern[]])?unitaryQ) ^ k_Integer := 
	unitaryPowerK[qop, k];
	
(* Composition of Unitary Operations *)
QuantumDiscreteOperation /: ((qop1:QuantumDiscreteOperation[specs1_Association, OptionsPattern[]])?unitaryQ)* ((qop2:QuantumDiscreteOperation[specs2_Association, OptionsPattern[]])?unitaryQ) := 
	unitaryComposition[qop1, qop2];
QuantumDiscreteOperation /: ((qop1:QuantumDiscreteOperation[specs1_Association, OptionsPattern[]])?unitaryQ)[((qop2:QuantumDiscreteOperation[specs2_Association, OptionsPattern[]])?unitaryQ)] := 
	unitaryComposition[qop1, qop2];
	
(* Symbolic Replacement for Quantum Operations *)	
QuantumDiscreteOperation /: ReplaceAll[qop:QuantumDiscreteOperation[specs_Association, OptionsPattern[]], subs_]:=
	operationReplaceAll[qop, subs];

(* Formatting Front End *)
QuantumDiscreteOperation /: MakeBoxes[qop:QuantumDiscreteOperation[specs_Association, OptionsPattern[]], fmt_]:=
	operationVisualize[qop, fmt];
		
(* Action of Quantum Operation on Quantum State *)
(qop:QuantumDiscreteOperation[specs_Association, OptionsPattern[]])[qstate_QuantumDiscreteState] := 
	qopAction[qop, qstate];



Clear[qopAction];

qopAction[qop_, qstate_] := With[{type = qop["OperationType"]},
	Which[
		type === "Unitary", unitaryAction[qop, qstate],
		type === "Oracle", oracleAction[qop, qstate],
		type === "BooleanFunction", boolAction[qop, qstate],
		type === "Projection", projectionAction[qop, qstate],
		type === "POVM", povmAction[qop, qstate],
		type === "Conditional", conditionAction[qop, qstate]]];


Clear[QuantumDiscreteOperationQ];
QuantumDiscreteOperationQ[expr_] := 
	If[Head[expr] === QuantumDiscreteOperation, True, False];


Clear[qopHelper];

(* Allowed types of quantum discrete operations *)
qopHelper["Unitary", spec_, dim_, pow_] := unitaryHelper[spec[[1]], spec[[2]], dim, pow];
qopHelper["Oracle", spec_, dim_, pow_] := oracleHelper[spec[[1]], spec[[2]], dim];
qopHelper["BooleanFunction", spec_, dim_, pow_] := boolFunHelper[spec[[1]], spec[[2]], dim];
qopHelper["Projection", spec_, dim_, pow_] := projectionHelper[spec[[1]], spec[[2]], dim];
qopHelper["POVM", spec_, dim_, pow_] := povmHelper[spec[[1]], spec[[2]], dim];
qopHelper["Conditional", spec_, dim_, pow_] := conditionHelper[spec];
qopHelper[type_, spec_, dim_, pow_] := (Message[QuantumDiscreteOperation::opType,type];$Failed);


ClearAll[opKeys, discreteOperationKeysQ];

opKeys := {
	"Operation",
	"OperationType",
	"Inputs",
	"QuditDimension",
	"Power"
	};

discreteOperationKeysQ[assoc_] := 
	If[Complement[Keys[assoc], opKeys] === {},
		True, False];


Clear[validDiscreteOperationQ];

validDiscreteOperationQ[qop_QuantumDiscreteOperation[assoc_Association, OptionsPattern[]]] := 
	If[discreteOperationKeysQ[qop[[1]]], True, False];
	
validDiscreteOperationQ[expr___] := False;


Clear[validDiscreteOperationListQ];

validDiscreteOperationListQ[ops_List] := Which[
	(*!AllTrue[ops, validDiscreteOperationQ],*) 
	!AllTrue[ops, QuantumDiscreteOperationQ],
	False,
	!AllTrue[ops[[All, 1, "QuditDimension"]], 
		# === ops[[1, 1, "QuditDimension"]]&],
	(Message[QuantumDiscreteOperation::listDims];False),
	True,
	True
	];
	
validDiscreteOperationListQ[expr___] := False;


(* ::Section:: *)
(*Symbolic Replacement for Quantum Operations*)


Clear[operationReplaceAll];

operationReplaceAll[qop_, subs_] := Module[
	{type, op, dim, assoc},
	type = qop["OperationType"];
	op = qop["Operation"];
	op = If[type === "Unitary" && sparseQ[op],
			SparseArray[ArrayRules[op]/.subs, Dimensions[op]],
			op/.subs];
		
	dim = qop["QuditDimension"];
	assoc = qop[[1]];
	assoc["Operation"] = op;
	QuantumDiscreteOperation[assoc, "QuditDimension" -> dim]
	];


(* ::Section:: *)
(*Operation Properties*)


Clear[extractQopProperty];

(* Ordered Input Quantum Objects *)
extractQopProperty[qop_, "Inputs"] :=  qop[[1, "Inputs"]];

(* Operation *)
extractQopProperty[qop_, "Operation"] := qop[[1, "Operation"]];

(* Operation Type *)
extractQopProperty[qop_, "OperationType"] := qop[[1, "OperationType"]];

(* Qudit Dimension *)
extractQopProperty[qop_, "QuditDimension"] := qop[[1, "QuditDimension"]];

(* Explicit Representation *)
extractQopProperty[qop_, "Representation"] := qopRepresentation[qop];

(* Hermitian *)
extractQopProperty[qop_, "HermitianQ"] := 
	HermitianMatrixQ @ qopRepresentation[qop];

(* Get inverse and conjugate transpose *)
extractQopProperty[qop_, "Inverse"] := qopInverse[qop];

extractQopProperty[qop_, "ConjugateTranspose"] := qopConjTranspose[qop];


extractQopProperty[expr___] := $Failed;


Clear[qopRepresentation];

qopRepresentation[qop_] := With[{type = qop[[1,"OperationType"]]},
	Which[
		type === "Unitary", unitaryRepresentation[qop],
		type === "Oracle", oracleRepresentation[qop],
		type === "BooleanFunction", boolRepresentation[qop],
		type === "Projection", projectionRepresentation[qop],
		type === "POVM", povmRepresentation[qop],
		type === "Conditional", conditionRepresentation[qop]]];


Clear[qopInverse];

qopInverse[qop_] := With[{type = qop[[1,"OperationType"]]},
	Which[
		type === "Unitary", unitaryInverse[qop],
		type === "Oracle", oracleInverse[qop],
		type === "BooleanFunction", boolInverse[qop],
		type === "Projection", projectionInverse[qop],
		type === "POVM", povmInverse[qop],
		type === "Conditional", conditionInverse[qop]]];


Clear[qopConjTranspose];

qopConjTranspose[qop_] := With[{type = qop[[1,"OperationType"]]},
	Which[
		type === "Unitary", unitaryConjugateTranspose[qop],
		type === "Oracle", oracleConjugateTranspose[qop],
		type === "BooleanFunction", boolConjugateTranspose[qop],
		type === "Projection", projectionConjugateTranspose[qop],
		type === "POVM", povmConjugateTranspose[qop],
		type === "Conditional", conditionConjugateTranspose[qop]]];


(* ::Section:: *)
(*Operation Front End Formatting*)


Clear[operationVisualize];

operationVisualize[qop_QuantumDiscreteOperation, fmt_] := 
	BoxForm`ArrangeSummaryBox[QuantumDiscreteOperation,
		qop,
		operationSpecVisual[qop],
		operationBaseVisual[qop], 
		operationExpandedVisual[qop], 
		fmt];


Clear[operationBaseVisual];
ClearAll[
	operationTypeVisualize, 
	operationDimVisualize,
	operationInputsVisualize
	];

operationBaseVisual[qop_] := Join[
	operationTypeVisualize[qop], 
	operationDimVisualize[qop],
	operationInputsVisualize[qop]
	];

operationTypeVisualize[qop_] := 
	{BoxForm`MakeSummaryItem[
		{"Operation Type: ",qop[[1,"OperationType"]]},
		 StandardForm]};

operationDimVisualize[qop_] := 
	{BoxForm`MakeSummaryItem[
		{"Qudit Dimension: ",qop[[1,"QuditDimension"]]},
		 StandardForm]};
		 
operationInputsVisualize[qop_] := If[
	qop[[1,"OperationType"]] != "Conditional",
	{BoxForm`MakeSummaryItem[
		{"Inputs: ",qop[[1,"Inputs"]]},
		 StandardForm]},
	{}];


Clear[operationSpecVisual];

operationSpecVisual[qop_] := With[{
	type = qop[[1,"OperationType"]]},
	Which[
		type === "Unitary", unitarySpecVisual[qop],
		type === "Oracle", oracleSpecVisual[qop],
		type === "BooleanFunction", boolSpecVisual[qop],
		type === "Projection", projectionSpecVisual[qop],
		type === "POVM", povmSpecVisual[qop],
		type === "Conditional", conditionSpecVisual[qop]]];


Clear[operationExpandedVisual];

operationExpandedVisual[qop_] := With[{
	type = qop[[1,"OperationType"]]},
	Which[
		type === "Unitary", unitaryExpandedVisual[qop],
		type === "Oracle", oracleExpandedVisual[qop],
		type === "BooleanFunction", boolExpandedVisual[qop],
		type === "Projection", projectionExpandedVisual[qop],
		type === "POVM", povmExpandedVisual[qop],
		type === "Conditional", conditionExpandedVisual[qop]]];


(* ::Subchapter:: *)
(*Unitary Operations*)


ClearAll[paulix, pauliy, pauliz, rootNot, hadamard, Rx, Ry, Rz, identity, fourier];
ClearAll[sigmap, sigmam];


sigmap := {{0, 1}, {0, 0}};
sigmam := {{0, 0}, {1, 0}};
pauliy := SparseArray[{{1,2} -> -I, {2,1} -> I}];
not := SparseArray[{{1,2} -> 1, {2,1} -> 1}];
rootNot := MatrixPower[not, 1/2];
hadamard := SparseArray[{{1,1} -> 1/Sqrt[2], {1,2} -> 1/Sqrt[2], {2,1} -> 1/Sqrt[2], {2,2} -> -1/Sqrt[2]}];

Rz[\[Theta]_] := SparseArray[{{1,1} -> 1, {2,2} -> Exp[I \[Theta]]}];
Rx[\[Theta]_] := {{Cos[\[Theta]/2], I*Sin[\[Theta]/2]}, {I * Sin[\[Theta]/2], Cos[\[Theta]/2]}};
Ry[\[Theta]_] := {{Cos[\[Theta]/2], Sin[\[Theta]/2]}, {- Sin[\[Theta]/2], Cos[\[Theta]/2]}};
paulix[d_, k_] := MatrixPower[SparseArray[({i_,j_}/;Mod[i-1,d,1]==j)->  1, {d,d}], k];
pauliz[d_, k_] := MatrixPower[SparseArray[{j_,j_}:> Exp[2 Pi I j/d], {d,d}], k];
identity[dim_] := SparseArray[{i_,i_}-> 1, {dim,dim}];
fourier[dim_, k_] := SparseArray[({i_,j_} :> Exp[ 2 Pi I Mod[(i-1)(j-1)*k,dim]/dim]/Sqrt[dim]), {dim,dim}];


ClearAll[sum, swap, rootSwap, control, cnot, cphase, toffoli, fredkin, deutsch, randUn];

sum[d_Integer] := SparseArray[{in_, out_} :> With[
	{i1 = IntegerDigits[in - 1, d, 2][[1]],
	 j1 = IntegerDigits[in - 1, d, 2][[2]],
	 i2 = IntegerDigits[out - 1, d, 2][[1]],
	 j2 = IntegerDigits[out - 1, d, 2][[2]]},
	 If[i1 == i2 && j2 == Mod[i1 + j1, d], 1, 0]], {d^2, d^2}];
	
swap[d_Integer] := SparseArray[({i_, j_}/; IntegerDigits[j - 1, d, 2] 
	== Reverse @ IntegerDigits[i - 1, d, 2]) -> 1, {d^2, d^2}];
	
rootSwap[d_Integer] := MatrixPower[swap[d], 1/2];

randUn[dim_Integer] := Module[
	{rr, rc},
	rr := RandomReal[NormalDistribution[0,1]];
	rc := rr + I rr;
	Orthogonalize[Table[rc, dim, dim]]
	];

control[dim_Integer, u_] := SparseArray[{i_, j_} :> Which[
	 IntegerDigits[i - 1, dim, 2][[1]] == 1 && IntegerDigits[j - 1, dim, 2][[1]] == 1,
	 u[[IntegerDigits[i - 1, dim, 2][[2]] + 1, IntegerDigits[j - 1, dim, 2][[2]] + 1]],
	 (i == j && IntegerDigits[i - 1, dim, 2][[1]] == 0), 1,
	 True, 0], {dim^2, dim^2}];
	 
cnot[dim_Integer] := SparseArray[{i_, j_} :> If[(IntegerDigits[j - 1, dim, 2] == 
	{IntegerDigits[i - 1, dim, 2][[1]], Mod[- Total @ IntegerDigits[i - 1, dim, 2], dim]}),
	1, 0], {dim^2, dim^2}];

cphase[d_Integer] := SparseArray[{i_, j_} :> Which[(i == j && IntegerDigits[i - 1, d, 2][[1]] == 0), 1,
	 i == j && IntegerDigits[j - 1, d, 2][[1]] > 0 && IntegerDigits[i - 1, d, 2][[1]] > 0, Exp[ 2 * Pi * I * (IntegerDigits[i - 1, d, 2][[2]])* (IntegerDigits[j - 1, d, 2][[2]]) /d],
	 True, 0], {d^2, d^2}];
	 
toffoli[nQu_Integer] := SparseArray[{i_, j_} :> 
	If[(i == j && i < 2^nQu - 1) || (i == 2^nQu - 1 && j == 2^nQu)|| (j == 2^nQu - 1 && i == 2^nQu), 1, 0], {2^nQu,2^nQu}];
	
fredkin = SparseArray[{{1,1} -> 1, {2,2} -> 1, {3,3} -> 1, {4,4} -> 1, {5,5} -> 1, {6,7} -> 1,
	{7,6} -> 1, {8,8} -> 1}, {8,8}];
	
deutsch[\[Theta]_] := SparseArray[{{1,1} -> 1, {2,2} -> 1, {3,3} -> 1, {4,4} -> 1, {5,5} -> 1, {6,6} -> 1,
	{7,7} -> I Cos[\[Theta]], {7,8} -> Sin[\[Theta]], {8,7} -> Sin[\[Theta]], {8,8} -> I Cos[\[Theta]]}, {8,8}];


Clear[unitaryExpandedVisual];
ClearAll[
	unitaryRepVisualize,
	hermQVisualize
	];
	
unitaryExpandedVisual[qop_] := Join[
	hermQVisualize[qop],
	unitaryRepVisualize[qop]
	];

unitaryRepVisualize[qop_] := With[{
	rep = unitaryRepresentation[qop]},
	If[
		Length @ rep <=  8,
		{BoxForm`MakeSummaryItem[
			{"Representation: ",
			 MatrixForm[Normal[rep]]},
			StandardForm]},
		{}]
	];
	
hermQVisualize[qop_] := With[{
	rep = unitaryRepresentation[qop]},
	{BoxForm`MakeSummaryItem[
		{"HermitianQ: ",
		 HermitianMatrixQ[rep]},
		 StandardForm]}];


Clear[unitaryPowerK];

unitaryPowerK[qop_QuantumDiscreteOperation, k_Integer] := Module[
	{assoc, dim, pow},
	assoc = qop[[1]];
	dim = assoc["QuditDimension"];
	pow = assoc["Power"];
	pow = Mod[pow * k, dim];
	assoc["Power"] = pow;
	QuantumDiscreteOperation[assoc, "QuditDimension" -> dim, "Power" -> pow]
	];


Clear[unitaryRepresentation];

unitaryRepresentation[qop_] := Module[
	{specs, op, ws, dim, pow, mat, perm, shape, size},
	specs = qop[[1]];
	op = specs["Operation"];
	ws = specs["Inputs"];
	dim = specs["QuditDimension"];
	pow = specs["Power"];
	If[unitarySingleQuditOperationQ[op["Base"]], 
		Return[unitarySingleQuditOperation[op, pow, dim]]]; 
	
	mat = unitaryMultiQuditOperation[op, pow, dim, Length @ ws];
	shape = ConstantArray[dim, 2 * Length @ ws];
	mat = ArrayReshape[mat, shape];
	size = dim ^ (Length @ ws);
	perm = InversePermutation[Ordering @ ws];
	perm = Join[perm, perm + Length @ perm];
	mat = Transpose[mat, perm];
	mat = ArrayReshape[Flatten @ mat, {size, size}];
	mat
	]


Clear[unitaryInverse];

unitaryInverse[qop_] := Module[
	{specs, op, invQ, dim, pow, base},
	specs = qop[[1]];
	op = qop["Operation"];
	
	base = op["Base"];
	
	invQ = op["InverseQ"];
	op["InverseQ"] = If[invQ, False, True];
	specs["Operation"] = op;
	dim = qop["QuditDimension"];
	pow = qop["Power"];
	
	QuantumDiscreteOperation[specs, "QuditDimension" -> dim, "Power" -> pow]
	];


Clear[unitaryConjugateTranspose];

unitaryConjugateTranspose[qop_] := Module[
	{specs, op, conjTranspQ, dim, pow},
	specs = qop[[1]];
	op = qop["Operation"];
	conjTranspQ = op["ConjugateTransposeQ"];
	op["ConjugateTransposeQ"] = If[conjTranspQ, False, True];
	specs["Operation"] = op;
	dim = qop["QuditDimension"];
	pow = qop["Power"];
	
	QuantumDiscreteOperation[specs, "QuditDimension" -> dim, "Power" -> pow]
	];


Clear[unitarySingleQuditOperationQ];

unitarySingleQuditOperationQ["SigmaX"|"SigmaY"|"SigmaZ"] := True;
unitarySingleQuditOperationQ["SigmaPlus"|"SigmaMinus"] := True;

unitarySingleQuditOperationQ[{"RotX"|"RotY"|"RotZ", \[Theta]_Number|\[Theta]_Real|\[Theta]_Rational|\[Theta]_Integer|\[Theta]_Symbol}] := True;
unitarySingleQuditOperationQ["Hadamard"|"NOT"|"RootNOT"|"S"|"T"] := True;

unitarySingleQuditOperationQ[expr___] := False;


Clear[unitarySingleQuditOperation];

unitarySingleQuditOperation[op_Association, k_Integer, d_Integer] := Module[
	{mat, invQ, conjTranspQ},
	mat = unitarySingleQuditOperationHelper[op["Base"], k, d];
	invQ = op["InverseQ"];
	conjTranspQ = op["ConjugateTransposeQ"];
	
	If[invQ, mat = Inverse[mat]];
	If[conjTranspQ, mat = ConjugateTranspose[mat]];
	mat
	];


Clear[unitarySingleQuditOperationHelper];

unitarySingleQuditOperationHelper["SigmaX", k_Integer, d_Integer] := paulix[d, k];
unitarySingleQuditOperationHelper["SigmaZ", k_Integer, d_Integer] := pauliz[d, k];

unitarySingleQuditOperationHelper[{"RotX", \[Theta]_Number|\[Theta]_Real|\[Theta]_Rational|\[Theta]_Integer|\[Theta]_Symbol}, 1, 2] := Rx[\[Theta]];
unitarySingleQuditOperationHelper[{"RotY", \[Theta]_Number|\[Theta]_Real|\[Theta]_Rational|\[Theta]_Integer|\[Theta]_Symbol}, 1, 2] := Ry[\[Theta]];
unitarySingleQuditOperationHelper[{"RotZ", \[Theta]_Number|\[Theta]_Real|\[Theta]_Rational|\[Theta]_Integer|\[Theta]_Symbol}, 1, 2] := Rz[\[Theta]];

unitarySingleQuditOperationHelper["SigmaY", 1, 2] := pauliy;
unitarySingleQuditOperationHelper["SigmaPlus", 1, 2] := sigmap;
unitarySingleQuditOperationHelper["SigmaMinus", 1, 2] := sigmam;

unitarySingleQuditOperationHelper["Hadamard", 1, 2] := hadamard;
unitarySingleQuditOperationHelper["NOT", 1, 2] := not;
unitarySingleQuditOperationHelper["RootNOT", 1, 2] := rootNot;
unitarySingleQuditOperationHelper["S", 1, 2] := Rz[Pi/2];
unitarySingleQuditOperationHelper["T", 1, 2] := Rz[Pi/4];


Clear[unitaryMultiQuditOperationQ];


unitaryMultiQuditOperationQ["CNOT"|"CPhase"] := True;
unitaryMultiQuditOperationQ["SUM"|"SWAP"|"RootSWAP"|"Fredkin"|"Toffoli"] := True;
unitaryMultiQuditOperationQ[{"Deutsch", \[Theta]_Number|\[Theta]_Real|\[Theta]_Rational|\[Theta]_Integer|\[Theta]_Symbol}] := True;

unitaryMultiQuditOperationQ["RandomUnitary"] := True;
unitaryMultiQuditOperationQ["Fourier"] := True;
unitaryMultiQuditOperationQ["Identity"] := True;

unitaryMultiQuditOperationQ[mat_?UnitaryMatrixQ] := True;
unitaryMultiQuditOperationQ[{"CU", u_?UnitaryMatrixQ}] := True;

unitaryMultiQuditOperationQ[expr___] := False;


Clear[unitaryMultiQuditOperation];

unitaryMultiQuditOperation[op_Association, k_Integer, d_Integer, n_Integer] := Module[
	{mat, invQ, conjTranspQ},
	mat = unitaryMultiQuditOperationHelper[op["Base"], k, d, n];
	invQ = op["InverseQ"];
	conjTranspQ = op["ConjugateTransposeQ"];
	
	If[invQ, mat = Inverse[mat]];
	If[conjTranspQ, mat = ConjugateTranspose[mat]];
	mat
	];


Clear[unitaryMultiQuditOperationHelper];

unitaryMultiQuditOperationHelper["RandomUnitary", 1, d_Integer, n_Integer] := randUn[d^n];
unitaryMultiQuditOperationHelper["Fourier", k_Integer, d_Integer, n_Integer] := fourier[d^n, k];
unitaryMultiQuditOperationHelper["Identity", 1, dim_Integer, nQu_Integer] := identity[dim^nQu];

unitaryMultiQuditOperationHelper["CNOT", 1, d_Integer, 2] := cnot[d];
unitaryMultiQuditOperationHelper["CPhase", 1, d_Integer, 2] := cphase[d];
unitaryMultiQuditOperationHelper["SUM", 1, d_Integer, 2] := sum[d];
unitaryMultiQuditOperationHelper["SWAP", 1, d_Integer, 2] := swap[d];
unitaryMultiQuditOperationHelper["RootSWAP", 1, d_Integer, 2] := rootSwap[d];
unitaryMultiQuditOperationHelper["Fredkin", 1, 2, 3] := fredkin;
unitaryMultiQuditOperationHelper["Toffoli", 1, 2, n_Integer] := toffoli[n];
unitaryMultiQuditOperationHelper[{"Deutsch", \[Theta]_Number|\[Theta]_Real|\[Theta]_Rational|\[Theta]_Integer|\[Theta]_Symbol}, 1, 2, 3] := deutsch[\[Theta]];

unitaryMultiQuditOperationHelper[{"CU", u_?UnitaryMatrixQ}, 1, d_Integer, 2] := control[d, u];
unitaryMultiQuditOperationHelper[mat_?UnitaryMatrixQ, k_Integer, d_Integer, n_Integer] := mat;


Clear[unitarySpecVisual];
unitarySpecVisual[qop_] := Module[
	{base, invQ, conjTranspQ, d, k, visual},
	base = unitaryBaseString[qop[[1,"Operation"]]];
	d = qop[[1,"QuditDimension"]];
	k = qop[[1, "Power"]];
	invQ = qop[[1, "Operation", "InverseQ"]];
	conjTranspQ =  qop[[1, "Operation", "ConjugateTransposeQ"]];
	
	visual = Which[d > 2 && k > 1,
		 ToString[Subsuperscript[base, d, k], StandardForm],
		 d > 2,
		 ToString[Subscript[base, d], StandardForm],
		 True,
		 base];
	If[ListQ[qop[[1, "Operation","Base"]]], 
		 Return[visual]];
		 
	visual = Which[
		(*ListQ[qop[[1, "Operation","Base"]]]&& MemberQ[{"RotX", "RotY", "RotZ"}, qop[[1, "Operation","Base",1]]], visual,*)
		invQ,ToString[Superscript[visual, "(-1)"], StandardForm],
		conjTranspQ, ToString[Superscript[visual, \[Dagger]], StandardForm],
		True, visual];
	
	visual
		 ]


Clear[unitaryHelper];

unitaryHelper[spec_, wires_, dim_, pow_] := Module[
	{ws, specFormatted, data},
	ws =  wiresFormat[wires];
	specFormatted = unitaryFormatHelper[spec, Length @ ws, dim];
	data = <|{"OperationType" -> "Unitary"}|>;
	data["Inputs"] = ws;
	data["QuditDimension"] = dim;
	data["Power"] = pow;
	data["Operation"] = 
		<|{"Base" -> specFormatted, "InverseQ" -> False, "ConjugateTransposeQ" -> False}|>;
	QuantumDiscreteOperation[data, "QuditDimension" -> dim, "Power" -> pow]
	];


Clear[unitaryFormatHelper];

unitaryFormatHelper["Identity"|"identity"|"id"|"Id", n_Integer, d_Integer] := "Identity";

unitaryFormatHelper["S-"|"SigmaMinus"|"sigmaMinus", n_Integer, d_Integer] := "SigmaMinus";
unitaryFormatHelper["S+"|"SigmaPlus"|"sigmaPlus", n_Integer, d_Integer] := "SigmaPlus";

unitaryFormatHelper["x"|"X"|"Sx"|"PauliX"|"SigmaX"|"pauliX"|"sigmaX", n_Integer, d_Integer] := "SigmaX";
unitaryFormatHelper["y"|"Y"|"Sy"|"PauliY"|"SigmaY"|"pauliY"|"sigmaY", n_Integer, d_Integer] := "SigmaY";
unitaryFormatHelper["z"|"Z"|"Sz"|"PauliZ"|"SigmaZ"|"pauliZ"|"sigmaZ", n_Integer, d_Integer] := "SigmaZ";

unitaryFormatHelper[{"Rx"|"RotX", \[Theta]_Number|\[Theta]_Real|\[Theta]_Rational|\[Theta]_Integer|\[Theta]_Symbol}, n_Integer, d_Integer] := {"RotX", \[Theta]};
unitaryFormatHelper[{"Ry"|"RotY", \[Theta]_Number|\[Theta]_Real|\[Theta]_Rational|\[Theta]_Integer|\[Theta]_Symbol}, n_Integer, d_Integer] := {"RotY", \[Theta]};
unitaryFormatHelper[{"Rz"|"RotZ", \[Theta]_Number|\[Theta]_Real|\[Theta]_Rational|\[Theta]_Integer|\[Theta]_Symbol}, n_Integer, d_Integer] := {"RotZ", \[Theta]};

unitaryFormatHelper["H"|"Hadamard"|"hadamard", n_Integer, d_Integer] := "Hadamard";
unitaryFormatHelper["Fourier"|"fourier"|"qft"|"QFT", n_Integer, d_Integer] := "Fourier";

unitaryFormatHelper["RootNOT"|"RootNot", n_Integer, d_Integer] := "RootNOT";
unitaryFormatHelper["NOT"|"Not"|"not", n_Integer, d_Integer] := "NOT";
unitaryFormatHelper["CNOT"|"Cnot"|"CNot"|"cnot", n_Integer, d_Integer] := "CNOT";
unitaryFormatHelper["CPhase"|"CPHASE"|"cphase"|"Cphase", n_Integer, d_Integer] := "CPHASE";
unitaryFormatHelper["SUM"|"sum"|"Sum", n_Integer, d_Integer] := "SUM";

unitaryFormatHelper["RootSWAP"|"RootSwap", n_Integer, d_Integer] := "RootSWAP";
unitaryFormatHelper["SWAP"|"swap"|"Swap", n_Integer, d_Integer] := "SWAP";

unitaryFormatHelper[{"Deutsch"|"deutsch", \[Theta]_Number|\[Theta]_Real|\[Theta]_Rational|\[Theta]_Integer|\[Theta]_Symbol}, n_Integer, d_Integer] := {"Deutsch", \[Theta]};
unitaryFormatHelper["Fredkin"|"fredkin"|"CSWAP"|"CSwap", n_Integer, d_Integer] := "Fredkin";
unitaryFormatHelper["Toffoli"|"toffoli", n_Integer, d_Integer] := "Toffoli";

unitaryFormatHelper["S", n_Integer, d_Integer] := "S";
unitaryFormatHelper["T", n_Integer, d_Integer] := "T";

unitaryFormatHelper["RandomUnitary", n_Integer, d_Integer] := randUn[d^n];
unitaryFormatHelper[u_?UnitaryMatrixQ, n_Integer, d_Integer] := u;
unitaryFormatHelper[{"CU"|"Cu"|"cu"|"Control"|"CONTROL"|"Controlled"|"CONTROLLED", u_}, n_Integer, d_Integer] := 
	If[UnitaryMatrixQ[u], 
		{"CU", u}, 
		(Message[QuantumDiscreteOperation::unitary,u];$Failed)];

unitaryFormatHelper[expr_, n_Integer, d_Integer] := (Message[QuantumDiscreteOperation::unitarySpec,expr];$Failed);


Clear[unitaryBaseString];

unitaryBaseString[assoc_] := 
	If[ListQ[assoc["Base"]] &&
		MemberQ[{"RotX", "RotY", "RotZ"}, assoc["Base"][[1]]],
		rotBaseString[assoc],
		unitaryBaseStringHelper[assoc["Base"]]];


Clear[rotBaseString];

rotBaseString[assoc_] := Module[
	{base, dir, angle, invQ, conjTranspQ, string},
	invQ = assoc["InverseQ"];
	conjTranspQ = assoc["ConjugateTransposeQ"];
	
	base = assoc["Base"];
	dir = Which[ 
		base[[1]] === "RotX", "x",
		base[[1]] === "RotY", "y",
		base[[1]] === "RotZ", "z"
		];
		
	angle = base[[2]];
	angle = If[Xor[invQ, conjTranspQ], -angle, angle];
	string = ToString[Subscript["R",dir], StandardForm]<>"("<>ToString[angle]<>")";
	string
];


Clear[unitaryBaseStringHelper];

unitaryBaseStringHelper["Identity"] := "\[ScriptCapitalI]\[ScriptD]";

unitaryBaseStringHelper["SigmaMinus"] := "\!\(\*SubscriptBox[\(\[Sigma]\), \(-\)]\)";
unitaryBaseStringHelper["SigmaPlus"] := "\!\(\*SubscriptBox[\(\[Sigma]\), \(+\)]\)";

unitaryBaseStringHelper["SigmaX"] := "\[ScriptCapitalX]";
unitaryBaseStringHelper["SigmaY"] := "\[ScriptCapitalY]";
unitaryBaseStringHelper["SigmaZ"] := "\[ScriptCapitalZ]";


unitaryBaseStringHelper["Hadamard"] := "\[ScriptCapitalH]";
unitaryBaseStringHelper["Fourier"] := "\[ScriptCapitalQ]\[ScriptCapitalF]\[ScriptCapitalT]";

unitaryBaseStringHelper["RootNOT"] := ToString[Sqrt["\[ScriptCapitalX]"], StandardForm];
unitaryBaseStringHelper["NOT"] := "\[ScriptCapitalX]";
unitaryBaseStringHelper["CNOT"] := ToString[Subscript["\[ScriptCapitalC]","\[ScriptCapitalX]"], StandardForm];
unitaryBaseStringHelper["CPHASE"] := ToString[Subscript["\[ScriptCapitalC]","\[ScriptCapitalZ]"], StandardForm];
unitaryBaseStringHelper["SUM"] := "\[ScriptCapitalS]\[ScriptU]\[ScriptM]";

unitaryBaseStringHelper["RootSWAP"] := ToString[Sqrt["\[ScriptCapitalS]\[ScriptCapitalW]\[ScriptCapitalA]\[ScriptCapitalP]"], StandardForm];
unitaryBaseStringHelper["SWAP"] := "\[ScriptCapitalS]\[ScriptCapitalW]\[ScriptCapitalA]\[ScriptCapitalP]";

unitaryBaseStringHelper[{"Deutsch", \[Theta]_Number|\[Theta]_Real|\[Theta]_Rational|\[Theta]_Integer|\[Theta]_Symbol}] := 
	"D("<>ToString[\[Theta]]<>")";
unitaryBaseStringHelper["Fredkin"] := "\[ScriptCapitalF]\[ScriptR]\[ScriptE]\[ScriptD]";
unitaryBaseStringHelper["Toffoli"] := "\[ScriptCapitalT]\[ScriptO]\[ScriptF]\[ScriptF]";

unitaryBaseStringHelper["S"] := "\[ScriptCapitalS]";
unitaryBaseStringHelper["T"] := "\[ScriptCapitalT]";

unitaryBaseStringHelper[u_?UnitaryMatrixQ] := "\[ScriptCapitalU]";
unitaryBaseStringHelper[{"CU", u_}] := ToString[Subscript["\[ScriptCapitalC]","\[ScriptCapitalU]"], StandardForm];


Clear[unitaryQ];

unitaryQ[qop_QuantumDiscreteOperation] := 
	If[qop[[1, "OperationType"]] === "Unitary", True, False];
	
unitaryQ[expr___] := False;


Clear[validUnitaryOperationQ];

validUnitaryOperationQ[qop_QuantumDiscreteOperation] := With[{
	baseOp = qop[[1,"Operation"]],
	pow = qop[[1,"Power"]],
	dim = qop[[1,"QuditDimension"]],
	arity = Length @ qop[[1,"Inputs"]]},
	
	unitaryOpQ[baseOp] &&
	unitaryArityQ[baseOp, arity] &&
	unitaryWiresQ[baseOp, dim] &&
	unitaryPowerQ[baseOp, pow, dim]];
	
validUnitaryOperationQ[expr___] := False;


Clear[unitaryOpQ];

unitaryOpQ["Identity"]:=True;
unitaryOpQ["RandomUnitary"]:=True;

unitaryOpQ["SigmaX"|"SigmaY"|"SigmaZ"]:=True;
unitaryOpQ["SigmaPlus"|"SigmaMinus"]:=True;

unitaryOpQ["Hadamard"]:=True;
unitaryOpQ["Fourier"]:=True;

unitaryOpQ["S"|"T"]:=True;
unitaryOpQ[{"RotX"|"RotY"|"RotZ",\[Theta]_Number|\[Theta]_Real|\[Theta]_Rational|\[Theta]_Integer|\[Theta]_Symbol}]:=True;

unitaryOpQ["NOT"|"RootNOT"]:=True;
unitaryOpQ["CNOT"|"CPHASE"|"SUM"]:=True;

unitaryOpQ["Toffoli"|"Fredkin"]:=True;
unitaryOpQ["SWAP"|"RootSWAP"]:=True;

unitaryOpQ[{"Deutsch", \[Theta]_Number|\[Theta]_Real|\[Theta]_Rational|\[Theta]_Integer|\[Theta]_Symbol}]:=True;

unitaryOpQ[{"CU", u_?UnitaryMatrixQ}]:=True;
unitaryOpQ[mat_?UnitaryMatrixQ] := True;

unitaryOpQ[expr___]:=(Message[QuantumDiscreteOperation::op,expr];False);


Clear[unitaryArityQ];

unitaryArityQ["Identity",1]:=True;
unitaryArityQ["SigmaX"|"SigmaY"|"SigmaZ"|"SigmaPlus"|"SigmaMinus",1]:=True;

unitaryArityQ["Hadamard"|"Fourier",1]:=True;

unitaryArityQ["S"|"T",1]:=True;
unitaryArityQ["NOT"|"RootNOT",1]:=True;
unitaryArityQ[{"RotX"|"RotY"|"RotZ", \[Theta]_Number|\[Theta]_Real|\[Theta]_Rational|\[Theta]_Integer|\[Theta]_Symbol},1]:=True;

unitaryArityQ["CNOT"|"SUM"|"CPHASE",2]:=True;
unitaryArityQ["SWAP"|"RootSWAP",2]:=True;

unitaryArityQ["Fredkin"|{"Deutsch", \[Theta]_Number|\[Theta]_Real|\[Theta]_Rational|\[Theta]_Integer|\[Theta]_Symbol}, 3]:=True;

unitaryArityQ["Toffoli", n_Integer]:=True;

unitaryArityQ[un_?UnitaryMatrixQ, n_Integer] := IntegerQ @ Power[Length @ un, 1/n];
unitaryArityQ[{"CU", u_?UnitaryMatrixQ}, 2] := True;

unitaryArityQ[expr__]:=(Message[QuantumDiscreteOperation::arity,expr];False);


Clear[unitaryWiresQ];

unitaryWiresQ["SigmaX"|"SigmaZ"|"Identity", d_Integer]:=True;
unitaryWiresQ["SigmaY"|"SigmaPlus"|"SigmaMinus", 2]:=True;

unitaryWiresQ["Hadamard", 2]:=True;
unitaryWiresQ["Fourier", d_Integer]:=True;
unitaryWiresQ["NOT"|"RootNOT", 2]:=True;
unitaryWiresQ[{"RotX"|"RotY"|"RotZ", \[Theta]_Number|\[Theta]_Real|\[Theta]_Rational|\[Theta]_Integer|\[Theta]_Symbol}, 2]:=True;
unitaryWiresQ["S"|"T", 2]:=True;

unitaryWiresQ["CNOT"|"CPHASE", d_Integer]:=True;
unitaryWiresQ["SUM", d_Integer]:=True;
unitaryWiresQ["RootSWAP", d_Integer]:=True;
unitaryWiresQ["SWAP", d_Integer]:=True;

unitaryWiresQ["Toffoli", 2]:=True;
unitaryWiresQ["Fredkin", 2]:=True;
unitaryWiresQ[{"Deutsch", \[Theta]_Number|\[Theta]_Real|\[Theta]_Rational|\[Theta]_Integer|\[Theta]_Symbol}, 2]:=True;

unitaryWiresQ[un_?UnitaryMatrixQ, d_Integer] := True;
unitaryWiresQ[{"CU", u_?UnitaryMatrixQ}, d_Integer]:=True;

unitaryWiresQ[expr___]:=(Message[QuantumDiscreteOperation::wires,expr];False);


Clear[unitaryPowerQ];

unitaryPowerQ["SigmaX"|"Fourier"|"SigmaZ", k_Integer, d_Integer]:= True;
unitaryPowerQ["Hadamard"|"SigmaY"|"SigmaPlus"|"SigmaMinus"|"Identity", 1, 2] := True;
unitaryPowerQ["S"|"T"|"NOT"|"RootNOT", 1, 2]:=True;
unitaryPowerQ[{"RotX"|"RotY"|"RotZ", \[Theta]_Number|\[Theta]_Real|\[Theta]_Rational|\[Theta]_Integer|\[Theta]_Symbol}, 1, 2]:=True;

unitaryPowerQ["CNOT"|"CPHASE"|"SUM", 1, d_Integer] := True;
unitaryPowerQ["SWAP"|"RootSWAP",1, d_Integer] := True;
unitaryPowerQ["Fredkin"|{"Deutsch", \[Theta]_Number|\[Theta]_Real|\[Theta]_Rational|\[Theta]_Integer|\[Theta]_Symbol}, 1, 2]:=True;
unitaryPowerQ["Toffoli", 1, 2]:=True;

unitaryPowerQ[un_?UnitaryMatrixQ, 1, d_Integer] := True;
unitaryPowerQ[{"CU", u_?UnitaryMatrixQ}, 1, d_Integer] := True;

unitaryPowerQ[expr__]:=(Message[QuantumDiscreteOperation::power,expr];False);


Clear[unitaryAction];

unitaryAction[qop_, qstate_] := Module[
	{opObjs, stateObjs, opDim, stateDim, pow, op, mat},
	opObjs = qop["Inputs"];
	stateObjs = qstate["QuantumObjects"];
	
	opDim = qop["QuditDimension"];
	stateDim = qstate["QuditDimension"];
	
	pow = qop["Power"];
	op = qop["Operation"];
	
	If[!AllTrue[opObjs, MemberQ[stateObjs, #]&],
		(Message[QuantumDiscreteOperation::opWiresNotInState,opObjs,qop,qstate];
		Return[$Failed])];
	If[opDim != stateDim,
		(Message[QuantumDiscreteOperation::opStateDims,opDim,stateDim];
		Return[$Failed])];
		
	mat = If[Length @ opObjs === 1,
		unitarySingleQuditMatrix[qop, opObjs[[1]], pow, opDim, qstate],
		unitaryMultiQuditMatrix[qop, opObjs, pow, opDim, qstate]
		];
	applyMatToState[mat, qstate]
	]


Clear[unitarySingleQuditMatrix];

unitarySingleQuditMatrix[qop_, obj_, pow_, dim_, qstate_] := Module[
	{stateObjs, opMats, pos},
	stateObjs = qstate["QuantumObjects"];
	opMats = ConstantArray[identity[dim], Length @ stateObjs];
	pos = Position[stateObjs, obj][[1,1]];
	opMats[[pos]] = qopRepresentation[qop];
	If[Length @ opMats == 1, opMats[[1]], KroneckerProduct @@ opMats]
	];


Clear[unitaryMultiQuditMatrix];

unitaryMultiQuditMatrix[qop_, objs_, pow_, dim_, qstate_] := Module[
	{stateWs, n, passiveWs, order, restOfRep, shape, opMat, tpLevels},
	stateWs = qstate["QuantumObjects"];
	n = Length @ stateWs;
	passiveWs = Complement[stateWs, objs];
	order = Join[objs, passiveWs];
	restOfRep = ConstantArray[identity[dim], Length @ passiveWs];
	shape = ConstantArray[dim, 2 * n];
	opMat = If[Length @ objs === n,
		qopRepresentation[qop],
		KroneckerProduct @@ Join[{qopRepresentation[qop]}, restOfRep]];
	opMat = ArrayReshape[opMat, shape];
	tpLevels = InversePermutation[Ordering @ order];
	tpLevels = Join[tpLevels, tpLevels + n];
	opMat = Transpose[opMat, tpLevels];
	opMat = ArrayReshape[Flatten @ opMat, {dim^n, dim^n}];
	opMat
	];


Clear[applyMatToState];

applyMatToState[mat_, state_] := Module[
	{newState},
	newState = state;
	If[vecQ[state], 
		newState[[1,"StateVector"]] = SparseArray[mat.state["StateVector"]], 
		newState[[1,"DensityMatrix"]] = SparseArray[mat.state["DensityMatrix"].ConjugateTranspose[mat]]
		];
	newState
	]


(* ::Subchapter:: *)
(*Oracles*)


Clear[oracleQ];

oracleQ[qop_QuantumDiscreteOperation] := 
	If[qop[[1, "OperationType"]] === "Oracle", True, False];
	
oracleQ[expr___] := False;


Clear[booleanFunctionQ];

booleanFunctionQ[qop_QuantumDiscreteOperation] := 
	If[qop[[1, "OperationType"]] === "BooleanFunction", True, False];
	
booleanFunctionQ[expr___] := False;


Clear[projectionQ];

projectionQ[qop_QuantumDiscreteOperation] := 
	If[qop[[1, "OperationType"]] === "Projection", True, False];
	
projectionQ[expr___] := False;


Clear[povmQ];

povmQ[qop_QuantumDiscreteOperation] := 
	If[qop[[1, "OperationType"]] === "POVM", True, False];
	
povmQ[expr___] := False;


Clear[conditionalQ];

conditionalQ[qop_QuantumDiscreteOperation] := 
	If[qop[[1, "OperationType"]] === "Conditional", True, False];
	
conditionalQ[expr___] := False;


(* ::Chapter:: *)
(*QuantumCircuit*)


Clear[QuantumCircuit];

(* Options *)
(*Options[QuantumCircuit] = {"Simplify" \[Rule] False};*)

(* Extract Properties of the Circuit *)
qcirc_QuantumCircuit[prop_String?StringQ] := extractQcircProperty[qcirc, prop];

(* Transform input into a normal form *)

(* From Operations *)
QuantumCircuit /: QuantumCircuit[ops_?validDiscreteOperationListQ] := 
	qcircFromOps[ops];
	
(* From Circuit *)
QuantumCircuit /: QuantumCircuit[circ_QuantumCircuit] := 
	qcircFromCircs[{circ}];

(* From Multiple Circuits *)
QuantumCircuit /: QuantumCircuit[circs_?circListQ] := 
	qcircFromCircs[circs];

QuantumCircuit /: circ1_QuantumCircuit + circ2_QuantumCircuit := 
	qcircFromCircs[{circ1, circ2}];
	
(* Circuit Exponentiation *)
QuantumCircuit /: circ_QuantumCircuit?validCircQ ^n_Integer :=
	qcircFromCircs[ConstantArray[circ, n]];
	
(* Symbolic Replacement in Circuit *)
QuantumCircuit /: ReplaceAll[circ_QuantumCircuit?validCircQ, subs_] :=
	qcircReplaceAll[circ, subs];
	
(* Formatting Front End *)
QuantumCircuit /: MakeBoxes[qcirc_QuantumCircuit?validCircQ, fmt_]:=
	circVisualize[qcirc, fmt];
		
(* Action of Quantum Circuit on Quantum State *)
(qcirc:qcirc_QuantumCircuit?validCircQ)[qstate_QuantumDiscreteState] := 
	qcircAction[qcirc, qstate];
	


Clear[qcircFromOps];

qcircFromOps[ops_, simp_] := Module[
	{assoc},
	assoc = <|{}|>;
	assoc["NumberOfOperations"] = Length @ ops;
	assoc["Operations"] = ops;
	assoc["Wires"] = Sort @ Union[Flatten @ Map[#["Inputs"]&, ops]];
	QuantumCircuit[assoc]
	];


Clear[qcircFromCircs];

qcircFromCircs[circs_] := Module[
	{assoc},
	assoc = <|{}|>;
	assoc["NumberOfOperations"] = Total @ Map[#["NumberOfOperations"]&, circs];
	assoc["Operations"] = Flatten @ Map[#["Operations"]&, circs];
	assoc["Wires"] = Sort @ Union @ Flatten @ Map[#["Wires"]&, circs];
	QuantumCircuit[assoc]
	];


Clear[validCircQ];

validCircQ[circ_] := 
	If[Head[circ] === QuantumCircuit &&
		AssociationQ @ circ[[1]],
		True, False]; 


Clear[circListQ];

circListQ[circs_] := 
	If[AllTrue[circs, Head[#] === QuantumCircuit &],
		True, False];


(* ::Section:: *)
(*Circuit Properties*)


Clear[extractQcircProperty];

(* Qudit Dimension *)
extractQcircProperty[qcirc_, "QuditDimension"] :=  qcirc[[1, "QuditDimension"]];

(* Circuit Wires *)
extractQcircProperty[qcirc_, "Wires"] :=  qcirc[[1, "Wires"]];

(* Operations *)
extractQcircProperty[qcirc_, "Operations"] := qcirc[[1, "Operations"]];

(* Number of Operations *)
extractQcircProperty[qcirc_, "NumberOfOperations"] := qcirc[[1, "NumberOfOperations"]];

extractQopProperty[expr___] := $Failed;


(* ::Section:: *)
(*Symbolic Replacement for Quantum Circuits*)


Clear[qcircReplaceAll];

qcircReplaceAll[circ_, subs_] := Module[
	{assoc, ops},
	assoc = circ[[1]];
	ops = circ["Operations"];
	ops = Map[operationReplaceAll[#, subs]&, ops];
	assoc["Operations"] = ops;
	QuantumCircuit[assoc]
	];


End[];
EndPackage[];


(* ::Input:: *)
(**)
