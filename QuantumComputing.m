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


QuantumFiniteDimensionalState::usage =
	"QuantumFiniteDimensionalState[{\!\(\*SubscriptBox[
StyleBox[\"spec\", \"TI\"], 
StyleBox[\"1\", \"TR\"]]\), \!\(\*SubscriptBox[
StyleBox[\"spec\", \"TI\"], \(2\)]\), \!\(\*
StyleBox[\"\[Ellipsis]\", \"TR\"]\),\*
StyleBox[\( \!\(\*
StyleBox[\" \", \"TR\"]\)\)]\!\(\*SubscriptBox[
StyleBox[\"spec\", \"TI\"], \(n\)]\)}] yields the finite dimensional "<>
	"quantum state resulting from taking the tensor product of the states specified by \!\(\*SubscriptBox[
StyleBox[\"spec\", \"TI\"], \(i\)]\).\n";
	
QuantumFiniteDimensionalStateQ::usage = 
	"QuantumFiniteDimensionalStateQ[\!\(\*
StyleBox[\"expr\", \"TI\"]\)] gives True if the head of \!\(\*
StyleBox[\"expr\", \"TI\"]\) is "<>
	"QuantumFiniteDimensionalState, and False otherwise.\n";
	
QuantumProduct::usage =
	"QuantumProduct[{\!\(\*SubscriptBox[
StyleBox[\"state\", \"TI\"], 
StyleBox[\"1\", \"TR\"]]\), \!\(\*SubscriptBox[
StyleBox[\"state\", \"TI\"], \(2\)]\), \!\(\*
StyleBox[\"\[Ellipsis]\", \"TR\"]\), \!\(\*
StyleBox[\" \", \"TR\"]\)\!\(\*SubscriptBox[
StyleBox[\"state\", \"TI\"], \(n\)]\)}] yields the finite dimensional "<>
	"quantum state resulting from taking the tensor product of \!\(\*SubscriptBox[
StyleBox[\"states\", \"TI\"], \(i\)]\).\n";

QuantumMixture::usage =
	"QuantumProduct[\!\(\*
StyleBox[\"states\", \"TI\"]\) -> \!\(\*
StyleBox[\"ws\", \"TI\"]\)] yields the finite dimensional "<>
	"quantum state generated from the statistical ensemble \!\(\*
StyleBox[\"states\", \"TI\"]\) weighted according to \!\(\*
StyleBox[\"ws\", \"TI\"]\).\n";
	
QuantumPartialTr::usage = 
	"QuantumPartialTr[\!\(\*
StyleBox[\"qstate\", \"TI\"]\), \!\(\*
StyleBox[\"qobs\", \"TI\"]\)\!\(\*
StyleBox[\"]\", \"TI\"]\)\!\(\*
StyleBox[\"  \", \"TI\"]\)yields the quantum state resulting from "<>
	"tracing out quantum objects \!\(\*
StyleBox[\"qobs\", \"TI\"]\) from quantum state \!\(\*
StyleBox[\"qstate\", \"TI\"]\)\n";
	
QuantumMatrixOperation::usage = 
	"QuantumMatrixOperation[\!\(\*
StyleBox[\"spec\", \"TI\"]\)] yields the quantum matrix "<>
	"operation with specification \!\(\*
StyleBox[\"spec\", \"TI\"]\)\!\(\*
StyleBox[\".\", \"TI\"]\)\!\(\*
StyleBox[\" \", \"TI\"]\)\n";

QuantumMatrixOperationQ::usage = 
	"QuantumMatrixOperationQ[\!\(\*
StyleBox[\"expr\", \"TI\"]\)] gives True if the head of \!\(\*
StyleBox[\"expr\", \"TI\"]\) is "<>
	"QuantumMatrixOperation, and False otherwise.\n";

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
StyleBox[\"2\", \"TR\"]]\)\!\(\*
StyleBox[\"]\", \"TI\"]\)\!\(\*
StyleBox[\" \", \"TI\"]\)yields the "<>
	"distance between quantum states \!\(\*SubscriptBox[
StyleBox[\"qstate\", \"TI\"], 
StyleBox[\"1\", \"TR\"]]\) and \!\(\*SubscriptBox[
StyleBox[\"qstate\", \"TI\"], 
StyleBox[\"2\", \"TR\"]]\).";
	
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
StyleBox[\"\[Ellipsis]\", \"TI\"]\) \n"<>
	
	"QuantumCircuit[\!\(\*
StyleBox[\"bf\", \"TI\"]\)] yields the reversible quantum circuit that "<>
	" generalizes the classical boolean function bf from bits to qubits.\n";

QuantumMeasurement::usage = 
	"QuantumMeasurement[\!\(\*
StyleBox[\"obs\", \"TI\"]\)\!\(\*
StyleBox[\"]\", \"TI\"]\)\!\(\*
StyleBox[\"  \", \"TI\"]\)yields the measurement operator whose eigenstates are given by projections onto "<>
	"the spectrum of observable operator \!\(\*
StyleBox[\"obs\", \"TI\"]\)\!\(\*
StyleBox[\".\", \"TI\"]\)\!\(\*
StyleBox[\" \", \"TI\"]\)\n"<>
	
	"QuantumMeasurement[\!\(\*
StyleBox[\"povm\", \"TI\"]\)\!\(\*
StyleBox[\"]\", \"TI\"]\)\!\(\*
StyleBox[\"  \", \"TI\"]\)yields the positive operator valued measurement \!\(\*
StyleBox[\"povm\", \"TI\"]\)\!\(\*
StyleBox[\".\", \"TI\"]\)\!\(\*
StyleBox[\" \", \"TI\"]\)\n";
	
QuantumMeasurementDistribution::usage = 
	"QuantumMeasurementDistribution[\!\(\*
StyleBox[\"qstate\", \"TI\"]\), \!\(\*
StyleBox[\"qmeas\", \"TI\"]\)\!\(\*
StyleBox[\" \", \"TI\"]\)\!\(\*
StyleBox[\"->\", \"TI\"]\)\!\(\*
StyleBox[\" \", \"TI\"]\)\!\(\*
StyleBox[\"ord\", \"TI\"]\)\!\(\*
StyleBox[\"]\", \"TI\"]\)\!\(\*
StyleBox[\"  \", \"TI\"]\)yields the symbolic distribution of eigenvalues resulting from "<>
	"the measurement of \!\(\*
StyleBox[\"qstate\", \"TI\"]\) according to QuantumMeasurement measurement operator \!\(\*
StyleBox[\"qmeas\", \"TI\"]\) applied with order \!\(\*
StyleBox[\"ord\", \"TI\"]\)\!\(\*
StyleBox[\".\", \"TI\"]\)\!\(\*
StyleBox[\" \", \"TI\"]\)\n";
	
QuantumEvaluate::usage = 
	"QuantumEvaluate[\!\(\*
StyleBox[\"qop\", \"TI\"]\) -> \!\(\*
StyleBox[\"order\", \"TI\"]\), \!\(\*
StyleBox[\"qstate\", \"TI\"]\)] evaluates the action of quantum matrix operation "<>
	"\!\(\*
StyleBox[\"qop\", \"TI\"]\)\!\(\*
StyleBox[\" \", \"TI\"]\)on quantum state \!\(\*
StyleBox[\"qstate\", \"TI\"]\)\!\(\*
StyleBox[\",\", \"TI\"]\) where \!\(\*
StyleBox[\"order\", \"TI\"]\) specifies which parts of \!\(\*
StyleBox[\"qstate\", \"TI\"]\) are acted upon.\n"<>
	
	"QuantumEvaluate[\!\(\*
StyleBox[\"qmeas\", \"TI\"]\) -> \!\(\*
StyleBox[\"order\", \"TI\"]\), \!\(\*
StyleBox[\"qstate\", \"TI\"]\)] evaluates the action of quantum measurement "<>
	"\!\(\*
StyleBox[\"qop\", \"TI\"]\)\!\(\*
StyleBox[\" \", \"TI\"]\)on quantum state \!\(\*
StyleBox[\"qstate\", \"TI\"]\)\!\(\*
StyleBox[\",\", \"TI\"]\) where \!\(\*
StyleBox[\"order\", \"TI\"]\) specifies which parts of \!\(\*
StyleBox[\"qstate\", \"TI\"]\) are acted upon.\n"<>
	
	"QuantumEvaluate[\!\(\*
StyleBox[\"qcirc\", \"TI\"]\) -> \!\(\*
StyleBox[\"order\", \"TI\"]\), \!\(\*
StyleBox[\"qstate\", \"TI\"]\)] evaluates quantum state \!\(\*
StyleBox[\"qstate\", \"TI\"]\)\!\(\*
StyleBox[\" \", \"TI\"]\)passing through quantum circuit "<>
	"\!\(\*
StyleBox[\"qcirc\", \"TI\"]\)\!\(\*
StyleBox[\",\", \"TI\"]\) where \!\(\*
StyleBox[\"order\", \"TI\"]\) specifies which parts of \!\(\*
StyleBox[\"qstate\", \"TI\"]\) are acted upon.\n";
	


(* ::Chapter:: *)
(*Function Definitions*)


Begin["`Private`"];


(* ::Section:: *)
(*QuantumFiniteDimensionalState*)


ClearAll[singleItemListQ, multiItemListQ];

singleItemListQ[lst_List]:= If[Length[lst] === 1, True, False];
singleItemListQ[expr___] := False;

multiItemListQ[lst_List]:= If[Length[lst] === 1, False, True];
multiItemListQ[expr___] := False;


Clear[ruleQ];
ruleQ[expr___] := If[Head[expr] === Rule, True, False];


Clear[stateArrayQ];
stateArrayQ[state_List] := 
	If[Depth[state] === 2 && 
		!StringQ[state[[1]]],
		True,
		False];
stateArrayQ[expr___] := False;


Clear[basisListQ];
basisListQ[arr_List] := 
	AllTrue[arr, # >= 0 && IntegerQ[#]&];
basisListQ[expr___] := False;







Clear[numQudits];
numQudits[{"BasisState", bases_}, d_Integer] := Length[bases];
numQudits[{"QuantumRegister", n_Integer, k_Integer}, d_Integer] := n;
numQudits["Plus", 2] := 1;
numQudits["Minus", 2] := 1;
numQudits["Left", 2] := 1;
numQudits["Right", 2] := 1;
numQudits["PsiPlus", 2] := 2;
numQudits["PhiPlus", 2] := 2;
numQudits["PsiMinus", 2] := 2;
numQudits["PhiMinus", 2] := 2;
numQudits[{"GHZ", n_Integer}, 2] := n;
numQudits[{"W", n_Integer}, 2] := n;
numQudits[{"RandomPure", n_Integer}, d_Integer] := n;
numQudits[{"PureUniformSuperposition", n_Integer}, d_Integer] := n;
numQudits[arr_?stateArrayQ, d_Integer] := Log[d, Length[arr]];






Clear[stateFormat];
stateFormat[{"BasisState", n_?Integer}] := {"BasisState", {n}};
stateFormat[{"QuantumRegister", n_Integer}] := {"QuantumRegister", n, 0};
stateFormat["+"|"plus"|"PLUS"] := "Plus";
stateFormat["-"|"minus"|"MINUS"] := "Minus";
stateFormat["l"|"L"|"left"|"LEFT"] := "Left";
stateFormat["r"|"R"|"right"|"RIGHT"] := "Right";
stateFormat["\[Psi]+"|"\[CapitalPsi]+"|"\!\(\*SubscriptBox[\(\[Psi]\), \(+\)]\)"|"\!\(\*SubscriptBox[\(\[CapitalPsi]\), \(+\)]\)"|"psiplus"|"PSIPLUS"] := "PsiPlus";
stateFormat["\[Phi]+"|"\[CapitalPhi]+"|"\!\(\*SubscriptBox[\(\[Phi]\), \(+\)]\)"|"\!\(\*SubscriptBox[\(\[CapitalPhi]\), \(+\)]\)"|"phiplus"|"PHIPLUS"] := "PhiPlus";
stateFormat["\[Psi]-"|"\[CapitalPsi]-"|"\!\(\*SubscriptBox[\(\[Psi]\), \(-\)]\)"|"\!\(\*SubscriptBox[\(\[CapitalPsi]\), \(-\)]\)"|"psiminus"|"PSIMINUS"] := "PsiMinus";
stateFormat["\[Phi]-"|"\[CapitalPhi]-"|"\!\(\*SubscriptBox[\(\[Phi]\), \(-\)]\)"|"\!\(\*SubscriptBox[\(\[CapitalPhi]\), \(-\)]\)"|"phiminus"|"PHIMINUS"] := "PhiMinus";
stateFormat["Ghz"|"ghz"|"GHZ"] := {"GHZ", 3};
stateFormat["PureUniformSuperposition"] := {"PureUniformSuperposition", 1};
stateFormat[{"Ghz", n_Integer}|{"ghz", n_Integer}] := {"GHZ", n};
stateFormat["w"|"W"] := {"W", 3};
stateFormat[{"w", n_Integer}|{"W", n_Integer}] := {"W", n};
stateFormat[expr___] := expr;


Clear[stateSpecFormat];
stateSpecFormat[spec_, dim_] := Module[
	{n, s},
	s = stateFormat[spec];
	n = numQudits[s, dim];
	{n, s}
	];






Clear[validStateSpecQHelper];
validStateSpecQHelper[{"BasisState", bases_?basisListQ}] := True;
validStateSpecQHelper[{"QuantumRegister", n_Integer, k_Integer}] := True;
validStateSpecQHelper[{"QuantumRegister", n_Integer}] := True;
validStateSpecQHelper["Plus"] := True;
validStateSpecQHelper["Minus"] := True;
validStateSpecQHelper["Left"] := True;
validStateSpecQHelper["Right"] := True;
validStateSpecQHelper["PsiPlus"] := True;
validStateSpecQHelper["PhiPlus"] := True;
validStateSpecQHelper["PsiMinus"] := True;
validStateSpecQHelper["PhiMinus"] := True;
validStateSpecQHelper[{"GHZ", n_Integer}] := True;
validStateSpecQHelper[{"W", n_Integer}] := True;
validStateSpecQHelper[{"RandomPure", n_Integer}] := True;
validStateSpecQHelper[{"PureUniformSuperposition", n_Integer}] := True;
validStateSpecQHelper[arr_?stateArrayQ] := True;
validStateSpecQHelper[expr___] := False;


Clear[validStateSpecQ];
validStateSpecQ[spec_] := 
	validStateSpecQHelper[
		stateFormat[spec]];








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
quantumRegister[num_, n_, d_] := 
	SparseArray[{{num+1} -> 1},{d^n}];






Clear[stateVec];
stateVec[n_, state_, dim_] := 
	FullSimplify @ Normalize @ SparseArray @ stateVecHelper[state, n, dim];


Clear[stateVecHelper];
(* Computational Basis States *)
stateVecHelper[{"BasisState", bases_}, n_Integer, d_Integer] := basisState[bases, d];
(* Quantum Registers *)
stateVecHelper[{"QuantumRegister", n_Integer, k_Integer}, n_Integer, d_Integer] :=
	quantumRegister[k, n, d];
(* +/-/L/R Basis States *)
stateVecHelper["Plus", 1, 2] := plus;
stateVecHelper["Minus", 1, 2] := minus;
stateVecHelper["Left", 1, 2] := left;
stateVecHelper["Right", 1, 2] := right;
(* Bell States *)
stateVecHelper["PsiPlus", 2, 2] := psiPlus;
stateVecHelper["PsiMinus", 2, 2] := psiMinus;
stateVecHelper["PhiPlus", 2, 2] := phiPlus;
stateVecHelper["PhiMinus", 2, 2] := phiMinus;
(* GHZ/W States *)
stateVecHelper[{"GHZ", n_Integer}, n_Integer, 2] := ghzState[n];
stateVecHelper[{"W", n_Integer}, n_Integer, 2] := wState[n];
(* Random Pure States *)
stateVecHelper[{"RandomPure", n_Integer}, n_Integer, d_Integer] := randPure[n, d];
(* Pure Uniform Superposition *)
stateVecHelper[{"PureUniformSuperposition", n_Integer}, n_Integer, d_Integer] := pureSuperpos[n, d];
(* States from Input Arrays *)
stateVecHelper[arr_?stateArrayQ, n_Integer, d_Integer] := FullSimplify @ Normalize @ arr;




ClearAll[vecToMat, matToVecHelper, matToVec];
vecToMat[vec_] := SparseArray[ConjugateTranspose[{vec}]].SparseArray[{vec}];
matToVecHelper[row_, pos_] := If[!(row[[pos]] ===  0), {row, pos}, False];
matToVec[mat_] := Module[
	{row, pos},
	{row, pos} = 
		SelectFirst[
			MapIndexed[matToVecHelper[#1, #2[[1]]]&, mat], !(False===#)&];
	row/Sqrt[row[[pos]]]
	];
	
Clear[l1Normalize];
l1Normalize[vals_List] := Map[#/(Total @ vals)&, vals];







(* Quantum State Queries *)
Clear[QuantumFiniteDimensionalStateQ];
QuantumFiniteDimensionalStateQ[s_] := 
	If[Head[s] === QuantumFiniteDimensionalState, True, False];
QuantumFiniteDimensionalStateQ[expr___] := False;

Clear[quantumFiniteDimensionalStateListQ];
quantumFiniteDimensionalStateListQ[states_List] := 
	If[AllTrue[states, QuantumFiniteDimensionalStateQ], True, False];
quantumFiniteDimensionalStateListQ[expr___] := False;

Clear[explicitQuantumStateQ];
explicitQuantumStateQ[state_QuantumFiniteDimensionalState] := 
	If[AssociationQ[state[[1]]], True, False];
explicitQuantumStateQ[expr___] := $Failed;

Clear[explicitQuantumStateListQ];
explicitQuantumStateListQ[states_List] := 
	(quantumFiniteDimensionalStateListQ[states] &&
	AllTrue[states, explicitQuantumStateQ] &&
	(!Length[states]===0));

ClearAll[pureStateQ, mixedStateQ];
pureStateQ[state_QuantumFiniteDimensionalState] := 
	Which[
		!AssociationQ[state[[1]]], $Failed,
		MemberQ[Keys[state[[1]]], "StateVector"], True,
		True, False];
mixedStateQ[state_QuantumFiniteDimensionalState] := 
	Which[
		!AssociationQ[state[[1]]], $Failed,
		MemberQ[Keys[state[[1]]], "DensityMatrix"], True,
		True, False];

ClearAll[qubitStateQ, singleQuditStateQ, singleQubitStateQ, singleQubitStateListQ];
qubitStateQ[state:QuantumFiniteDimensionalState[assoc_Association]] := 
	If[assoc["QuditDimension"] === 2, True, False];
singleQuditStateQ[state:QuantumFiniteDimensionalState[assoc_Association]] := 
	If[assoc["NumberOfQudits"] === 1, True, False];
singleQubitStateQ[state:QuantumFiniteDimensionalState[assoc_Association]] := 
	qubitStateQ[state] && singleQuditStateQ[state];
singleQubitStateListQ[states_List] := 
	(quantumFiniteDimensionalStateListQ[states] && 
	AllTrue[states, singleQubitStateQ]);
















(* Extract Properties of the State *)

ClearAll[getStateVector, getDensityMatrix];
getStateVector[qs_] := If[pureStateQ[qs], qs[[1, "StateVector"]], 
		(Message[QuantumFiniteDimensionalState::mixture,qs];$Failed)];
getDensityMatrix[qs_] := If[mixedStateQ[qs], qs[[1, "DensityMatrix"]], 
		vecToMat[qs[[1, "StateVector"]]]];

ClearAll[getBasisLabels, getBasisStates];
getBasisLabels[qs_] := Range[qs[[1,"NumberOfQudits"]]] /. List -> Ket;
getBasisStates[qs_] := Module[
	{dim, n, kets},
	dim = qs[[1,"QuditDimension"]];
	n = qs[[1,"NumberOfQudits"]];
	kets = Tuples[Range[0, dim-1], n];
	kets = Map[ToString, kets, {2}];
	kets = Map[StringJoin, kets];
	Map[Ket, kets]
	];

ClearAll[getVonNeumannEntropy, getPurity];
getVonNeumannEntropy[qs_] := Module[
	{dm, eigs, entropy},
	If[pureStateQ[qs], Return[0]];
	dm = qs[[1,"DensityMatrix"]];
	eigs = Select[Eigenvalues[dm],#>0&];
	entropy = - Total@Map[# Log[2,#]&,eigs];
	entropy
	];
getPurity[qs_] := If[pureStateQ[qs], 
	1, With[{dm = qs[[1,"DensityMatrix"]]}, Abs[Tr[dm.dm]]]];

ClearAll[blochPureSphericalCoords, blochPureCartesianCoords];
ClearAll[blochMixedSphericalCoords, blochMixedCartesianCoords];
ClearAll[getBlochSphericalCoords, getBlochCartesianCoords];
ClearAll[u, v, w, r, \[Theta], \[Phi]];

blochPureSphericalCoords[qs_] := 
	With[{
		\[Alpha] = qs[[1,"StateVector"]][[1]],
		\[Beta] = qs[[1,"StateVector"]][[2]]
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
	dm = Normal @ qs[[1,"DensityMatrix"]];
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
	
getBlochSphericalCoords[qs_?singleQubitStateQ] := With[{
	vals = If[pureStateQ[qs], 
		blochPureSphericalCoords[qs],
		blochMixedSphericalCoords[qs]]},
	Thread[{"r", "\[Theta]", "\[Phi]"} -> vals]];
getBlochSphericalCoords[expr___] := $Failed;
	
getBlochCartesianCoords[qs_?singleQubitStateQ] := With[{
	vals = If[pureStateQ[qs], 
		blochPureCartesianCoords[qs],
		blochMixedCartesianCoords[qs]]},
	Thread[{"u", "v", "w"} -> vals]];
getBlochCartesianCoords[expr___] := $Failed;









(* Bloch Sphere Visualization *)
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

blochQubit[qs_QuantumFiniteDimensionalState, OptionsPattern[]] := With[{
	coords = getBlochCartesianCoords[qs][[All, 2]],
	color = OptionValue["Color"]},
	Graphics3D[{
		color, 
		Arrowheads[.03], 
		Arrow[Tube[{{0,0,0},coords},.01],
		{0, -.01}]},
		Boxed->False,
		PlotRange->{{-1.7, 1.7},{-1.7, 1.7},{-1.7, 1.7}}
		]
	];


Clear[colorListQ];
colorListQ[colors_List] := If[AllTrue[colors, ColorQ], True, False];
colorListQ[expr___] := False;


Clear[QuantumBlochPlot];

Options[QuantumBlochPlot] = {PlotStyle -> True};

QuantumBlochPlot[qstate_?singleQubitStateQ, OptionsPattern[]] := 
	QuantumBlochPlot[{qstate}, PlotStyle -> OptionValue[PlotStyle]];
QuantumBlochPlot[qstates_?singleQubitStateListQ, OptionsPattern[]] := Module[
	{cols, sLength, cLength},
	cols = OptionValue[PlotStyle];
	sLength = Length[qstates];
	cLength = Length[cols];
	cols = Which[
		cols === True,
		RandomColor[sLength],
		cLength === sLength,
		cols,
		cLength > sLength,
		cols[[1;;sLength]],
		True,
		Join[cols, RandomColor[sLength - cLength]]];
	Show[Join[{greatCircles, referenceStates}, 
			MapThread[blochQubit[#1, "Color" -> #2]&, {qstates, cols}]]]
		];
		
		
(* Phase and Amplitude Plotting *)

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
	vec = qs[[1,"StateVector"]];
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
	mat = qs[[1,"DensityMatrix"]];
	
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
	
Clear[phaseAndAmplitudePlot];
phaseAndAmplitudePlot[qs4_] := If[pureStateQ[qs4], plotVec[qs4], plotMat[qs4]];
	
	
	


Clear[extractQStateProperty];

(* Basics *)
extractQStateProperty[qs_, "QuditDimension"] := qs[[1, "QuditDimension"]];
extractQStateProperty[qs_, "NumberOfQudits"] := qs[[1, "NumberOfQudits"]];
extractQStateProperty[qs_, "StateVector"] := getStateVector[qs];	
extractQStateProperty[qs_, "DensityMatrix"] := getDensityMatrix[qs];

(* Basis *)
extractQStateProperty[qs_, "BasisStates"] := getBasisStates[qs];

(* Numerical Properties *)
extractQStateProperty[qs_, "VonNeumannEntropy"] := getVonNeumannEntropy[qs];
extractQStateProperty[qs_, "Purity"] := getPurity[qs];

(* State Queries *)
extractQStateProperty[qs_, "PureStateQ"] := pureStateQ[qs];
extractQStateProperty[qs_, "MixedStateQ"] := mixedStateQ[qs];

(* Bloch Sphere *)
extractQStateProperty[qs_, "BlochCartesianCoordinates"] := getBlochCartesianCoords[qs];
extractQStateProperty[qs_, "BlochSphericalCoordinates"] := getBlochSphericalCoords[qs];
extractQStateProperty[qs_, "BlochPlot"] := QuantumBlochPlot[qs];

(* Visualizations *)
extractQStateProperty[qs_, "Plot"] := phaseAndAmplitudePlot[qs];

extractQStateProperty[expr___] := $Failed;



(* Quantum Finite Dimensional States *)
ClearAll[QuantumFiniteDimensionalState, iQState];
QuantumFiniteDimensionalState::mixture="Mixed state `1` cannot be represented by a state vector";

(* Options *)
Options[QuantumFiniteDimensionalState] = {"QuditDimension" -> 2};

(* Transform input into a normal form *)
iQState[{numQudits_, state_}, dim_] := Module[
	{assoc},
	assoc = <|{
		"QuditDimension" -> dim, 
		"NumberOfQudits" -> numQudits,
		"StateVector" -> stateVec[numQudits, state, dim]}|>;
	
	QuantumFiniteDimensionalState[assoc]
	];


QuantumFiniteDimensionalState /: qstate:QuantumFiniteDimensionalState[spec_?singleItemListQ, OptionsPattern[]] := 
	If[validStateSpecQ[spec[[1]]],
		iQState[stateSpecFormat[spec[[1]], OptionValue["QuditDimension"]], OptionValue["QuditDimension"]], 
		HoldForm[qstate]];
QuantumFiniteDimensionalState /: QuantumFiniteDimensionalState[specs_?multiItemListQ, OptionsPattern[]] := 
	QuantumProduct[Map[QuantumFiniteDimensionalState[{#}, "QuditDimension" -> OptionValue["QuditDimension"]]&, specs]];

QuantumFiniteDimensionalState /: qstate_QuantumFiniteDimensionalState[prop_String?StringQ] := 
	extractQStateProperty[qstate, prop];


(* Front End Fornatting of State *)
QuantumFiniteDimensionalState/: MakeBoxes[qs3:QuantumFiniteDimensionalState[specs_Association, OptionsPattern[]], fmt_]:=
	stateVisualize[qs3, fmt];


Clear[stateVisualize];

stateVisualize[qs2_QuantumFiniteDimensionalState, fmt_] := 
	BoxForm`ArrangeSummaryBox[QuantumFiniteDimensionalState,
		qs2,
		stateArrayVisualize[qs2],
		stateBaseVisual[qs2], 
		stateExpandedVisual[qs2], 
		fmt];
		
		
Clear[stateArrayVisualize];
stateArrayVisualize[qs_] := Which[
	pureStateQ[qs] && Dimensions[qs[[1,"StateVector"]]][[1]] < 8, 
	MatrixForm @ Normal @ qs[[1,"StateVector"]],
	pureStateQ[qs], 
	qs[[1,"StateVector"]],
	True,
	qs[[1,"DensityMatrix"]]]
	
	
Clear[stateBaseVisual];
ClearAll[stateObjsVisualize, stateDimVisualize];
stateBaseVisual[qs1_] := 
	Join[stateObjsVisualize[qs1], stateDimVisualize[qs1]];
	
stateObjsVisualize[qs_] := 
	{BoxForm`MakeSummaryItem[
		{"Number Of Qudits: ",qs[[1,"NumberOfQudits"]]},
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

stateExpandedVisual[qs5_] := Join[
	purityVisualize[qs5],
	entropyVisualize[qs5],
	blochVisualize[qs5]
	];
	
purityVisualize[qs_] := 
	{BoxForm`MakeSummaryItem[
		{"Quantum Purity: ", getPurity[qs]},
		 StandardForm]};

entropyVisualize[qs_] := 
	{BoxForm`MakeSummaryItem[
		{"Von Neumann Entropy: ", getVonNeumannEntropy[qs]},
		 StandardForm]};

blochVisualize[qs_] := 
	If[singleQubitStateQ[qs],
		{BoxForm`MakeSummaryItem[
			{"Bloch Coordinates: ", getBlochSphericalCoords[qs]},
			 StandardForm]},
		{}
	];


(* Product States *)

ClearAll[tensorPure, tensorMixed];

tensorPure[subStates_, dim_] := Module[
	{rep, n, assoc},
	rep = Flatten[KroneckerProduct @@ subStates[[All, 1, "StateVector"]]];
	n = Total @ Map[#["NumberOfQudits"]&, subStates];
	assoc = <|{
		"StateVector" -> rep,
		"QuditDimension" -> dim,
		"NumberOfQudits" -> n
		}|>;
	QuantumFiniteDimensionalState[assoc]
	];
	
tensorMixed[subStates_, dim_] := Module[
	{rep, n, mats, assoc},
	mats = Map[If[pureStateQ[#], 
		vecToMat @ #[[1,"StateVector"]],
		#[[1, "DensityMatrix"]]]&, subStates];
	rep = KroneckerProduct @@ mats;
	n = Total @ Map[#["NumberOfQudits"]&, subStates];
	assoc = <|{
		"DensityMatrix" -> SparseArray[rep],
		"QuditDimension" -> dim,
		"NumberOfQudits" -> n
		}|>;
	QuantumFiniteDimensionalState[assoc]
	];

Clear[validStateProductListQ];
validStateProductListQ[qstates_List] := 
	If[
		quantumFiniteDimensionalStateListQ[qstates] &&
		AllTrue[qstates, AssociationQ[#[[1]]]&],
		True,
		False];
		
ClearAll[QuantumProduct, iQProduct];
QuantumProduct[states1_?validStateProductListQ] := iQProduct[states1];
iQProduct[subStates_] := Module[
	{dim},
	dim = subStates[[1]]["QuditDimension"];
	If[!AllTrue[subStates, #["QuditDimension"] === dim &],
		(Message[QuantumFiniteDimensionalState::prodDims,subStates];Return[$Failed])];
	If[AllTrue[subStates, pureStateQ[#]&], 
		tensorPure[subStates, dim],
		tensorMixed[subStates, dim]]
	];


(* Mixtures *)
ClearAll[QuantumMixture, iQMixture];
QuantumMixture[states2_ -> weights_] := 
	If[ListQ[weights] && explicitQuantumStateListQ[states2],
		iQMixture[states2, weights], HoldForm[QuantumMixture[states2 -> weights]]];
	
iQMixture[states3_, weights_] := Module[
	{probs, mixture, reps, dm, pureQ, hasStateVec, hasDM},
	If[Length[states3] === 1, Return[states3[[1]]]];
	If[!AllTrue[states3, #["NumberOfQudits"] === states3[[1]]["NumberOfQudits"]&],
		Return[$Failed]];
		
	probs = l1Normalize[weights];
	mixture = states3[[1, 1]];
	reps = Map[If[pureStateQ[#], 
		vecToMat @ #[[1,"StateVector"]],
		#[[1, "DensityMatrix"]]]&, states3];
		
	dm = Total @ MapThread[Times,{reps, probs}];
	
	mixture = If[MemberQ[Keys[mixture], "StateVector"],
		Delete[mixture, "StateVector"],
		Delete[mixture, "DensityMatrix"]];
	
	pureQ = If[Tr[dm.dm] === 1, True, False];
	If[pureQ, 
		mixture["StateVector"] = matToVec[dm],
		mixture["DensityMatrix"] = dm];	
	QuantumFiniteDimensionalState[mixture]
	];
iQMixture[expr___] := $Failed;



(* Tracing out subsystems *)
ClearAll[QuantumPartialTr, iQPartialTr];
QuantumPartialTr[s_, qsOut_] := 
	If[ListQ[qsOut] && QuantumFiniteDimensionalStateQ[s],
		iQPartialTr[s, qsOut], HoldForm[QuantumPartialTr[s, qsOut]]];

iQPartialTr[s_, {}] := s;
iQPartialTr[s_, qsOut_] := Module[
	{newState, dm, dim, wContract, n, tpLevels, assoc},
	dm = s["DensityMatrix"];
	dim = s["QuditDimension"];
	wContract = qsOut[[1]];
	n = s["NumberOfQudits"];
	dm = ArrayReshape[dm, ConstantArray[dim, 2 * n]];
	tpLevels = InversePermutation[Join[Range[n], Range[n] + n]];
	dm = Transpose[dm, tpLevels];
	dm = TensorContract[dm, {{wContract, wContract + n}}];
	tpLevels = InversePermutation@Join[Range[n - 1], Range[n - 1] + (n - 1)];
	dm = Transpose[dm, tpLevels];
	dm = ArrayReshape[Flatten @ dm, {dim^(n-1), dim^(n-1)}];
	assoc = <|{
		"NumberOfQudits" -> n - Length[qsOut],
		"QuditDimension" -> dim
		}|>;
	If[Tr[dm.dm] == 1,
		assoc["StateVector"] = matToVec[dm],
		assoc["DensityMatrix"] = dm];
	newState = QuantumFiniteDimensionalState[assoc];
	iQPartialTr[newState, Drop[qsOut, 1]]
	];
iQPartialTr[expr___] := $Failed;


(* ::Section:: *)
(*QuantumMatrixOperation*)


ClearAll[paulix, pauliy, pauliz, rootNot, hadamard, Rx, Ry, Rz, fourier, id];
ClearAll[sigmap, sigmam];

id[d_] := IdentityMatrix[d];
sigmap := {{0, 1}, {0, 0}};
sigmam := {{0, 0}, {1, 0}};
pauliy := SparseArray[{{1,2} -> -I, {2,1} -> I}];
not := SparseArray[{{1,2} -> 1, {2,1} -> 1}];
rootNot := MatrixPower[not, 1/2];
hadamard := SparseArray[{{1,1} -> 1/Sqrt[2], {1,2} -> 1/Sqrt[2], {2,1} -> 1/Sqrt[2], {2,2} -> -1/Sqrt[2]}];

Rz[\[Theta]1_] := SparseArray[{{1,1} -> 1, {2,2} -> Exp[I \[Theta]1]}];
Rx[\[Theta]2_] := {{Cos[\[Theta]2/2], I*Sin[\[Theta]2/2]}, {I * Sin[\[Theta]2/2], Cos[\[Theta]2/2]}};
Ry[\[Theta]3_] := {{Cos[\[Theta]3/2], Sin[\[Theta]3/2]}, {- Sin[\[Theta]3/2], Cos[\[Theta]3/2]}};
paulix[d1_, k1_] := MatrixPower[SparseArray[({i_,j_}/;Mod[i-1,d1,1]==j)->  1, {d1,d1}], k1];
pauliz[d_, k_] := MatrixPower[SparseArray[{j_,j_}:> Exp[2 Pi I j/d], {d,d}], k];
fourier[dim_, k_] := SparseArray[({i_,j_} :> Exp[ 2 Pi I Mod[(i-1)(j-1)*k,dim]/dim]/Sqrt[dim]), {dim,dim}];

ClearAll[sum, swap, rootSwap, control, cnot, cphase, toffoli, fredkin, deutsch, randUn];

sum[d_Integer] := SparseArray[{in_, out_} :> With[
	{i1 = IntegerDigits[in - 1, d, 2][[1]],
	 j1 = IntegerDigits[in - 1, d, 2][[2]],
	 i2 = IntegerDigits[out - 1, d, 2][[1]],
	 j2 = IntegerDigits[out - 1, d, 2][[2]]},
	 If[i1 == i2 && j2 == Mod[i1 + j1, d], 1, 0]], {d^2, d^2}];
	
swap[d1_Integer] := SparseArray[({i_, j_}/; IntegerDigits[j - 1, d1, 2] 
	== Reverse @ IntegerDigits[i - 1, d1, 2]) -> 1, {d1^2, d1^2}];
	
rootSwap[dim_Integer] := MatrixPower[swap[dim], 1/2];

randUn[dim_Integer, n_Integer] := Module[
	{rr, rc},
	rr := RandomReal[NormalDistribution[0,1]];
	rc := rr + I rr;
	Orthogonalize[Table[rc, dim^n, dim^n]]
	];

control[dim_Integer, u_] := SparseArray[{i_, j_} :> Which[
	 IntegerDigits[i - 1, dim, 2][[1]] == 1 && IntegerDigits[j - 1, dim, 2][[1]] == 1,
	 u[[IntegerDigits[i - 1, dim, 2][[2]] + 1, IntegerDigits[j - 1, dim, 2][[2]] + 1]],
	 (i == j && IntegerDigits[i - 1, dim, 2][[1]] == 0), 1,
	 True, 0], {dim^2, dim^2}];
	 
cnot[d2_Integer] := SparseArray[{i_, j_} :> If[(IntegerDigits[j - 1, d2, 2] == 
	{IntegerDigits[i - 1, d2, 2][[1]], Mod[- Total @ IntegerDigits[i - 1, d2, 2], d2]}),
	1, 0], {d2^2, d2^2}];

cphase[d_Integer] := SparseArray[{i_, j_} :> Which[(i == j && IntegerDigits[i - 1, d, 2][[1]] == 0), 1,
	 i == j && IntegerDigits[j - 1, d, 2][[1]] > 0 && IntegerDigits[i - 1, d, 2][[1]] > 0, Exp[ 2 * Pi * I * (IntegerDigits[i - 1, d, 2][[2]])* (IntegerDigits[j - 1, d, 2][[2]]) /d],
	 True, 0], {d^2, d^2}];
	 
toffoli[n_Integer] := SparseArray[{i_, j_} :> 
	If[(i == j && i < 2^n- 1) || (i == 2^n - 1 && j == 2^n)|| (j == 2^n - 1 && i == 2^n), 1, 0], {2^n,2^n}];
	
fredkin = SparseArray[{{1,1} -> 1, {2,2} -> 1, {3,3} -> 1, {4,4} -> 1, {5,5} -> 1, {6,7} -> 1,
	{7,6} -> 1, {8,8} -> 1}, {8,8}];
	
deutsch[\[Theta]_] := SparseArray[{{1,1} -> 1, {2,2} -> 1, {3,3} -> 1, {4,4} -> 1, {5,5} -> 1, {6,6} -> 1,
	{7,7} -> I Cos[\[Theta]], {7,8} -> Sin[\[Theta]], {8,7} -> Sin[\[Theta]], {8,8} -> I Cos[\[Theta]]}, {8,8}];
	
	

Clear[opArrayQ];
opArrayQ[op_List] := 
	If[Depth[op] === 3 && AllTrue[op, ListQ],
	True,
	False];
opArrayQ[expr___] := False;



Clear[opSpecFormat];
opSpecFormat["S-"|"SigmaMinus"|"sigmaMinus"] := "SigmaMinus";
opSpecFormat["S+"|"SigmaPlus"|"sigmaPlus"] := "SigmaPlus";
opSpecFormat["x"|"X"|"Sx"|"PauliX"|"SigmaX"|"pauliX"|"sigmaX"] := "SigmaX";
opSpecFormat["y"|"Y"|"Sy"|"PauliY"|"SigmaY"|"pauliY"|"sigmaY"] := "SigmaY";
opSpecFormat["z"|"Z"|"Sz"|"PauliZ"|"SigmaZ"|"pauliZ"|"sigmaZ"] := "SigmaZ";
opSpecFormat[{"Rx"|"RotX", \[Theta]_}] := {"RotX", \[Theta]};
opSpecFormat[{"Ry"|"RotY", \[Theta]_}] := {"RotY", \[Theta]};
opSpecFormat[{"Rz"|"RotZ", \[Theta]_}] := {"RotZ", \[Theta]};
opSpecFormat["H"|"Hadamard"|"hadamard"] := "Hadamard";
opSpecFormat["Fourier"|"fourier"|"qft"|"QFT"] := {"Fourier", 1};
opSpecFormat[{"fourier", n_Integer}] := {"Fourier", n};
opSpecFormat[{"QFT", n_Integer}] := {"Fourier", n};
opSpecFormat[{"qft", n_Integer}] := {"Fourier", n};
opSpecFormat["RootNOT"|"RootNot"] := "RootNOT";
opSpecFormat["NOT"|"Not"|"not"] := "NOT";
opSpecFormat["CNOT"|"Cnot"|"CNot"|"cnot"] := "CNOT";
opSpecFormat["CPhase"|"CPHASE"|"cphase"|"Cphase"] := "CPHASE";
opSpecFormat["SUM"|"sum"|"Sum"] := "SUM";
opSpecFormat["RootSWAP"|"RootSwap"] := "RootSWAP";
opSpecFormat["SWAP"|"swap"|"Swap"] := "SWAP";
opSpecFormat[{"Deutsch"|"deutsch", \[Theta]_}] := {"Deutsch", \[Theta]};
opSpecFormat["Fredkin"|"fredkin"|"CSWAP"|"CSwap"] := "Fredkin";
opSpecFormat["Toffoli"|"toffoli"] := {"Toffoli", 3};
opSpecFormat[{"toffoli", n_Integer}] := {"Toffoli", n};
opSpecFormat[{"Toffoli", n_Integer}] := {"Toffoli", n};
opSpecFormat["S"] := "S";
opSpecFormat["T"] := "T";
opSpecFormat["Tdag"] := "Tdag";
opSpecFormat["RandomUnitary"] := {"RandomUnitary", 2};
opSpecFormat[expr___] := expr;


Clear[validAngleQ];
validAngleQ[\[Theta]_Number|\[Theta]_Real|\[Theta]_Rational|\[Theta]_Integer|\[Theta]_Symbol|\[Theta]_Times] := True;
validAngleQ[expr___] := False;



Clear[validOpSpecQ];
validOpSpecQ[{"RandomUnitary", n_Integer}] := True;
validOpSpecQ["SigmaMinus"] := True;
validOpSpecQ["SigmaPlus"] := True;
validOpSpecQ["SigmaX"] := True;
validOpSpecQ["SigmaY"] := True;
validOpSpecQ["SigmaZ"] := True;
validOpSpecQ[{"RotX", \[Theta]_?validAngleQ}] := True;
validOpSpecQ[{"RotY", \[Theta]_?validAngleQ}] := True;
validOpSpecQ[{"RotZ", \[Theta]_?validAngleQ}] := True;
validOpSpecQ["Hadamard"] := True;
validOpSpecQ["Fourier"] := True;
validOpSpecQ[{"Fourier", n_Integer}] := True;
validOpSpecQ["RootNOT"] := True;
validOpSpecQ["NOT"] := True;
validOpSpecQ["CNOT"] := True;
validOpSpecQ["CPHASE"] := True;
validOpSpecQ["SUM"] := True;
validOpSpecQ["RootSWAP"] := True;
validOpSpecQ["SWAP"] := True;
validOpSpecQ[{"Deutsch", \[Theta]_?validAngleQ}] := True;
validOpSpecQ["Fredkin"] := True;
validOpSpecQ["Toffoli"] := True;
validOpSpecQ[{"Toffoli", n_Integer}] := True;
validOpSpecQ["S"] := True;
validOpSpecQ["T"] := True;
validOpSpecQ["Tdag"] := True;
validOpSpecQ[arr_?opArrayQ] := True;
validOpSpecQ[expr___] := False;




Clear[inGateSetQ];
inGateSetQ["SigmaX", 2, 1] := True;
inGateSetQ["SigmaY", 2, 1] := True;
inGateSetQ["SigmaZ", 2, 1] := True;
inGateSetQ["Hadamard", 2, 1] := True;
inGateSetQ["CNOT", 2, 1] := True;
inGateSetQ["T", 2, 1] := True;
inGateSetQ["Tdag", 2, 1] := True;
inGateSetQ[expr___] := False;




Clear[opMat];
opMat[{"RandomUnitary", n_Integer}, d_, pow_] := MatrixPower[randUn[d, n], pow];
opMat["SigmaX", d_, pow_] := paulix[d, pow];
opMat["SigmaZ", d_, pow_] := pauliz[d, pow];
opMat["SigmaY", d_, pow_] := MatrixPower[pauliy, pow];
opMat["SigmaPlus", d_, pow_] := MatrixPower[sigmap, pow];
opMat["SigmaMinus", d_, pow_] := MatrixPower[sigmap, pow];
opMat[{"RotX", \[Theta]_}, d_, pow_] := Rx[\[Theta] * pow];
opMat[{"RotY", \[Theta]_}, d_, pow_] := Ry[\[Theta] * pow];
opMat[{"RotZ", \[Theta]_}, d_, pow_] := Rz[\[Theta] * pow];
opMat["Hadamard", d_, pow_] := MatrixPower[hadamard, pow];
opMat["NOT", d_, pow_] := MatrixPower[hadamard, pow];
opMat["RootNOT", d_, pow_] := MatrixPower[not, pow/2];
opMat["S", d_, pow_] := Rz[pow * Pi/2];
opMat["T", d_, pow_] := Rz[pow * Pi/4];
opMat["Tdag", d_, pow_] := Rz[- pow * Pi/4];
opMat[{"Fourier", n_}, d_, pow_] := fourier[d^n, pow];
opMat["CNOT", d_, pow_] := MatrixPower[cnot[d], pow];
opMat["CPHASE", d_, pow_] := MatrixPower[cphase[d], pow];
opMat["SUM", d_, pow_] := MatrixPower[sum[d], pow];
opMat["SWAP", d_, pow_] := MatrixPower[swap[d], pow];
opMat["RootSWAP", d_, pow_] := MatrixPower[swap[d], pow/2];
opMat["Fredkin", d_, pow_] := MatrixPower[fredkin, pow];
opMat[{"Toffoli", n_}, d_, pow_] := MatrixPower[toffoli[n], pow];
opMat[{"Deutsch", \[Theta]_}, d_, pow_] := MatrixPower[deutsch[\[Theta]], pow];
opMat[expr_, d_, pow_] := MatrixPower[expr, pow];


Clear[opHasNameQ];
opHasNameQ[qop_] := MemberQ[Keys[qop[[1]]], "OperationName"];

Clear[notAssociationQ];
notAssociationQ[expr___] := If[AssociationQ[expr], False, True];


Clear[extractQMatProperty];
extractQMatProperty[qop_, "QuditDimension"] := qop[[1, "QuditDimension"]];
extractQMatProperty[qop_, "Arity"] := qop[[1, "Arity"]];
extractQMatProperty[qop_, "MatrixRepresentation"] := qop[[1, "Matrix"]];	
extractQMatProperty[qop_, "InGateSetQ"] := qop[[1, "InGateSetQ"]];
extractQMatProperty[qop_?opHasNameQ, "OperationName"] := qop[[1, "OperationName"]];
extractQMatProperty[qop_, "HermitianQ"] := HermitianMatrixQ[qop[[1, "Matrix"]]];
extractQMatProperty[qop_, "UnitaryQ"] := UnitaryMatrixQ[qop[[1, "Matrix"]]];
extractQMatProperty[expr___] := $Failed;


Clear[QuantumMatrixOperationQ];
QuantumMatrixOperationQ[expr___] := 
	If[Head[expr] === QuantumMatrixOperation, True, False];





		
Clear[opSpecVisualize];
opSpecVisualize[op_] := 
	If[opHasNameQ[op], op[[1]]["OperationName"], op[[1]]["Matrix"]];

ClearAll[opBaseVisualize, opDimVisualize, opArityVisualize];
opBaseVisualize[qop_] := Join[
	opDimVisualize[qop],
	opArityVisualize[qop]
	];
	
opDimVisualize[qop_] := 
	{BoxForm`MakeSummaryItem[
		{"Qudit Dimension: ",qop[[1,"QuditDimension"]]},
		 StandardForm]};

opArityVisualize[qop_] := 
	{BoxForm`MakeSummaryItem[
		{"Arity: ",qop[[1,"Arity"]]},
		 StandardForm]};
		 
ClearAll[opExpandedVisualize, hermQVisualize, unitaryQVisualize];
opExpandedVisualize[qop_] := Join[
	hermQVisualize[qop],
	unitaryQVisualize[qop]
	];
	
hermQVisualize[qop_] := With[
	{mat = qop[[1, "Matrix"]]},
	{BoxForm`MakeSummaryItem[
		{"Hermitian: ", If[HermitianMatrixQ[mat], 
		"True", "False"]},
		 StandardForm]}];
		 
unitaryQVisualize[qop_] := With[
	{mat = qop[[1, "Matrix"]]},
	{BoxForm`MakeSummaryItem[
		{"Unitary: ", If[UnitaryMatrixQ[mat], 
		"True", "False"]},
		 StandardForm]}];
		 
Clear[opMatVisualize];
opMatVisualize[qop2_QuantumMatrixOperation, fmt_] := 
	BoxForm`ArrangeSummaryBox[QuantumMatrixOperation,
		qop2,
		opSpecVisualize[qop2],
		opBaseVisualize[qop2], 
		opExpandedVisualize[qop2], 
		fmt];
		

	ClearAll[QuantumMatrixOperation, iQMatrix];

(* Options *)
Options[QuantumMatrixOperation] = {"QuditDimension" -> 2, "Power" -> 1};

(* Transform input into a normal form *)
iQMatrix[opSpec_, dim_, pow_] := Module[
	{mat, assoc, gateOpBool},
	mat = opMat[opSpec, dim, pow];
	assoc = <|{
		"QuditDimension" -> dim, 
		"Matrix" -> mat,
		"Arity" -> Log[dim, Length[mat]]}|>;
	gateOpBool = inGateSetQ[opSpec, dim, pow];
	assoc["InGateSetQ"] = gateOpBool;
	If[!stateArrayQ[opSpec], assoc["OperationName"] = 
		Which[
			StringQ[opSpec], 
			opSpec,
			ListQ[opSpec] && StringQ[opSpec[[1]]],
			opSpec[[1]]]
			];
	QuantumMatrixOperation[assoc]
	];

QuantumMatrixOperation /: qop:QuantumMatrixOperation[spec_?notAssociationQ, OptionsPattern[]]:=
	If[validOpSpecQ[opSpecFormat[spec]], 
		iQMatrix[opSpecFormat[spec], OptionValue["QuditDimension"], OptionValue["Power"]],
		HoldForm[qop]];
	
		(* Front End Fornatting of Mat Op *)
QuantumMatrixOperation/: MakeBoxes[qmat1:QuantumMatrixOperation[specs_Association, OptionsPattern[]], fmt_]:=
	opMatVisualize[qmat1, fmt];
		
QuantumMatrixOperation /: qmat_QuantumMatrixOperation[prop_String?StringQ] := 
	extractQMatProperty[qmat, prop];
	
Clear[explicitQuantumMatrixOpQ];
explicitQuantumMatrixOpQ[op_QuantumMatrixOperation] := 
	If[Length[op] === 1 && AssociationQ[op[[1]]], True, False];
explicitQuantumMatrixOpQ[expr___] := False;


(* ::Section:: *)
(*QuantumMeasurement*)


ClearAll[validObsSpecQ, validPOVMSpecQ];

validObsSpecQ["SigmaX"|"SigmaY"|"SigmaZ"|"Hadamard"] := True;
validObsSpecQ[expr_] := HermitianMatrixQ[expr];

validPOVMSpecQ[expr_] := With[{n = Length[expr[[1]]]},
	If[(Simplify @ Total @ Map[ConjugateTranspose[#].#&, expr]) === IdentityMatrix[n],
		True, False]];

Clear[inMeasSetQ];
inMeasSetQ["SigmaX"|"SigmaZ", 2] := True;
inMeasSetQ[expr___] := False;

Clear[obsMatrix];
obsMatrix["SigmaX", d_] := paulix[d, 1];
obsMatrix["SigmaY", 2] := pauliy;
obsMatrix["SigmaZ", d_] := pauliz[d, 1];
obsMatrix["Hadamard", 2] := hadamard;
obsMatrix[expr_] := expr;

Clear[hasEigsQ];
hasEigsQ[qmeas_] := If[MemberQ[Keys[qmeas[[1]]], "Eigenvalues"], True, False];

Clear[hasMsQ];
hasMsQ[qmeas_] := If[MemberQ[Keys[qmeas[[1]]], "Ms"], True, False];

Clear[extractQMeasProperty];
extractQMeasProperty[qmeas_, "QuditDimension"] := qmeas[[1, "QuditDimension"]];
extractQMeasProperty[qmeas_, "Arity"] := qmeas[[1, "Arity"]];
extractQMeasProperty[qmeas_, "MeasurementType"] := qmeas[[1, "MeasurementType"]];	
extractQMeasProperty[qmeas_, "InMeasurementSetQ"] := qmeas[[1, "InMeasurementSetQ"]];
extractQMeasProperty[qmeas_?opHasNameQ, "OperationName"] := qmeas[[1, "OperationName"]];
extractQMeasProperty[qmeas_?hasEigsQ, "Eigenvalues"] := qmeas[[1, "Eigenvalues"]];
extractQMeasProperty[qmeas_?hasEigsQ, "Eigenvectors"] := qmeas[[1, "Eigenvectors"]];
extractQMeasProperty[qmeas_?hasMsQ, "Ms"] := qmeas[[1, "Ms"]];
extractQMeasProperty[expr___] := $Failed;


Clear[QuantumMeasurementQ];
QuantumMeasurementQ[expr___] := 
	If[Head[expr] === QuantumMeasurement, True, False];



Clear[measSpecVisualize];
measSpecVisualize[meas1_] := 
	If[MemberQ[Keys[meas1[[1]]], "OperationName"], 
		meas1[[1]]["OperationName"],
		meas1[[1]]["MeasurementType"]
		];

ClearAll[measBaseVisualize, measDimVisualize, measArityVisualize];
measBaseVisualize[meas4_] := Join[
	measDimVisualize[meas4],
	measArityVisualize[meas4]
	];
	
measDimVisualize[meas_] := 
	{BoxForm`MakeSummaryItem[
		{"Qudit Dimension: ",meas[[1,"QuditDimension"]]},
		 StandardForm]};

measArityVisualize[meas_] := 
	{BoxForm`MakeSummaryItem[
		{"Arity: ",meas[[1,"Arity"]]},
		 StandardForm]};
		 
Clear[MeasVisualize];
MeasVisualize[qmeas2_QuantumMeasurement, fmt_] := 
	BoxForm`ArrangeSummaryBox[QuantumMeasurement,
		qmeas2,
		measSpecVisualize[qmeas2],
		measBaseVisualize[qmeas2], 
		{},
		fmt];
		 

		
		
		

ClearAll[QuantumMeasurement, iQObservable, iQPOVM];

iQObservable[obs_, d_] := Module[
	{mat, assoc, measOpBool, eigVals, eigVecs},
	mat = obsMatrix[obs, d];
	assoc = <|{
		"QuditDimension" -> d, 
		"MeasurementType" -> "Projection",
		"Matrix" -> mat,
		"Arity" -> Log[d, Length[mat]]}|>;
	measOpBool = inMeasSetQ[obs, d];
	assoc["InMeasurementSetQ"] = measOpBool;
	If[StringQ[obs], assoc["OperationName"] = obs];
	If[ListQ[obs] && StringQ[obs[[1]]],
		assoc["OperationName"] = obs[[1]]];
	
	{eigVals, eigVecs} = Eigensystem[mat];
	assoc["Eigenvalues"] = eigVals;
	assoc["Eigenvectors"] = Normalize /@ eigVecs;
	QuantumMeasurement[assoc]
	];
	
iQPOVM[povm_, d_] := Module[
	{assoc},
	assoc = <|{
		"QuditDimension" -> d, 
		"MeasurementType" -> "POVM",
		"InMeasurementSetQ" -> False,
		"Ms" -> povm,
		"Arity" -> Log[d, Length[povm[[1]]]]}|>;
	QuantumMeasurement[assoc]
	];

Options[QuantumMeasurement] = {"QuditDimension" -> 2};
QuantumMeasurement /: qmeas:QuantumMeasurement["Observable" -> obs_, OptionsPattern[]] :=
	If[validObsSpecQ[obs], 
		iQObservable[obs, OptionValue["QuditDimension"]],
		HoldForm[qmeas]];
QuantumMeasurement /: qmeas:QuantumMeasurement["POVM" -> povm_, OptionsPattern[]] :=
	If[validPOVMSpecQ[povm], 
		iQPOVM[povm, OptionValue["QuditDimension"]],
		HoldForm[qmeas]];
		
		(* Front End Fornatting of Measurement *)
QuantumMeasurement/: MakeBoxes[qmeas1:QuantumMeasurement[specs_Association, OptionsPattern[]], fmt_]:=
	MeasVisualize[qmeas1, fmt];
		
QuantumMeasurement /: qmeas_QuantumMeasurement[prop_String?StringQ] := 
	extractQMeasProperty[qmeas, prop];
	
Clear[explicitQuantumMeasurementQ];
explicitQuantumMeasurementQ[meas_QuantumMeasurement] := 
	If[Length[meas] === 1 && AssociationQ[meas[[1]]], True, False];
explicitQuantumMeasurementQ[expr___] := False;





	
	
	
Clear[singleQuditProjector];
singleQuditProjector[qmeas_, qs_, pos_] := Module[
	{projs, dim, reps},
	projs = eigToProjector /@ qmeas["Eigenvectors"];
	dim = qs["QuditDimension"];
	If[qs["NumberOfQudits"] === 1, 
		Return[projs]];
	reps = ConstantArray[id[dim], qs["NumberOfQudits"]];
	reps = Map[ReplacePart[reps, pos -> #]&, projs];
	reps = Map[KroneckerProduct @@ # &, reps];
	reps
];

Clear[multiQuditProjector];
multiQuditProjector[qmeas_, qs_, pos_] := Module[
	{projs, qObjs, arity, n, dim, passiveObjs, order, reps, passiveRep, shape, tpLevels},
	projs = eigToProjector /@ qmeas["Eigenvectors"];
	n = qs["NumberOfQudits"];
	dim = qs["QuditDimension"];
	
	qObjs = Range[n];
	arity = qmeas["Arity"];

	passiveObjs = Complement[qObjs, pos];
	order = Join[pos, passiveObjs];
	passiveRep = ConstantArray[id[dim], Length @ passiveObjs];
	shape = ConstantArray[dim, 2 * n];
	reps = If[Length @ passiveObjs === 0,
		projs,
		Map[KroneckerProduct @@ Join[{#}, passiveRep]&, projs]];	
	reps = Map[ArrayReshape[#, shape]&, reps];
	
	tpLevels = InversePermutation[Ordering @ order];
	tpLevels = Join[tpLevels, tpLevels + n];
	reps = Map[Transpose[#, tpLevels]&, reps];
	reps = Map[ArrayReshape[Flatten @ #, {dim^n, dim^n}]&, reps];
	reps
	];
	
Clear[singleQuditPOVM];
singleQuditPOVM[qmeas_, qs_, pos_] := Module[
	{povm, dim, reps},
	povm = Map[ConjugateTranspose[#].#&, qmeas["Ms"]];
	dim = qs["QuditDimension"];
	If[qs["NumberOfQudits"] === 1, 
		Return[povm]];
		
	reps = ConstantArray[id[dim], qs["NumberOfQudits"]];
	reps = Map[ReplacePart[reps, pos -> #]&, povm];
	reps = Map[KroneckerProduct @@ # &, reps];
	reps
];

Clear[multiQuditPOVM];
multiQuditPOVM[qmeas_, qs_, pos_] := Module[
	{povm, qObjs, arity, n, dim, passiveObjs, order, reps, passiveRep, shape, tpLevels},
	povm = Map[ConjugateTranspose[#].#&, qmeas["Ms"]];
	n = qs["NumberOfQudits"];
	dim = qs["QuditDimension"];
	
	qObjs = Range[n];
	arity = qmeas["Arity"];

	passiveObjs = Complement[qObjs, pos];
	order = Join[pos, passiveObjs];
	passiveRep = ConstantArray[id[dim], Length @ passiveObjs];
	shape = ConstantArray[dim, 2 * n];
	reps = If[Length @ passiveObjs === 0,
		povm,
		Map[KroneckerProduct @@ Join[{#}, passiveRep]&, povm]];	
	reps = Map[ArrayReshape[#, shape]&, reps];
	
	tpLevels = InversePermutation[Ordering @ order];
	tpLevels = Join[tpLevels, tpLevels + n];
	reps = Map[Transpose[#, tpLevels]&, reps];
	reps = Map[ArrayReshape[Flatten @ #, {dim^n, dim^n}]&, reps];
	reps
	];
	
Clear[qMeasProbsAndVals];
qMeasProbsAndVals[qs_, qm_, ord_] := Module[
	{dm, mType, arity, probs, vals},
	dm = qs["DensityMatrix"];
	mType = qm["MeasurementType"];
	vals = If[mType === "Projection",
		qm["Eigenvalues"],
		Range[Length[qm["Ms"]]]];
	arity = qm["Arity"];
	probs = Which[
		arity === 1 && mType === "Projection",
		singleQuditProjector[qm, qs, ord[[1]]],
		arity === 1 && mType === "POVM",
		singleQuditPOVM[qm, qs, ord[[1]]],
		mType === "Projection",
		multiQuditProjector[qm, qs, ord],
		True,
		multiQuditPOVM[qm, qs, ord]];
	probs = Map[Tr[#.dm]&, probs];
	{probs, vals}
	];
	

Clear[QuantumMeasurementDistribution];
QuantumMeasurementDistribution[ord_ -> qm_QuantumMeasurement?explicitQuantumMeasurementQ, 
	qs_QuantumFiniteDimensionalState?explicitQuantumStateQ] := With[
	{probVals = qMeasProbsAndVals[qs, qm, ord]},
	EmpiricalDistribution[probVals[[1]] -> probVals[[2]]]];


(* ::Subsection:: *)
(*QuantumEntangledObjectsQ*)


Clear[QuantumEntangledObjectsQ];

QuantumEntangledObjectsQ::inputs="Arguments to QuantumEntangledObjectsQ must be "<>
	"a pure quantum finiteDimensional state qstate and a bipartition of the quantum objects "<>
	"in qstate";
QuantumEntangledObjectsQ::qstate="Argument `1` at position 1 is not a quantum finiteDimensional state";
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
	
	
Clear[validEntangledStateQ];

validEntangledStateQ[qs_QuantumFiniteDimensionalState] := If[
	qs["PureStateQ"],
	True,
	(Message[QuantumEntangledObjectsQ::pureState,qs];False)];

validEntangledStateQ[expr___] := 
	(Message[QuantumEntangledObjectsQ::qstate,expr];False);
	

Clear[validBipartitionQ];

validBipartitionQ[qs_, split_List] := Which[
	Length[split] != 2,
	(Message[QuantumEntangledObjectsQ::splitLength,split];False),
	!AllTrue[split[[1]],MemberQ[Range[qs["NumberOfQudits"]], #]&],
	(Message[QuantumEntangledObjectsQ::qObjs,qs,split[[1]],1,split];False),
	!AllTrue[split[[2]],MemberQ[Range[qs["NumberOfQudits"]], #]&],
	(Message[QuantumEntangledObjectsQ::qObjs,qs,split[[2]],2,split];False),
	!DuplicateFreeQ[Flatten @ split],
	(Message[QuantumEntangledObjectsQ::unique,split];False),
	!Complement[Range[qs["NumberOfQudits"]], Flatten@split] === {},
	(Message[QuantumEntangledObjectsQ::bipartition,qs,split];False),
	True,
	True
	];
	
validBipartitionQ[qs_, split_] := 
	(Message[QuantumEntangledObjectsQ::splitList,split];False);


(* ::Subsection:: *)
(*QuantumStateDistance *)


Clear[QuantumStateDistance];

QuantumStateDistance::inputs="Inputs to QuantumStateDistance must be two "<>
	"quantum finiteDimensional states";
QuantumStateDistance::distMeasure="Method `1` is not a valid option for "<>
	"calculating the distance between two states";
QuantumStateDistance::qubit="Argument `1` must be a single qubit state for "<>
	" Euclidean distance method";
QuantumStateDistance::dims="Arguments `1` at position 1 and `2` at position 2 "<>
	"must have the same qudit dimension";
QuantumStateDistance::numObjs="Arguments `1` at position 1 and `2` at position 2 "<>
	"must describe the same number of quantum objects";

(* Options *)
Options[QuantumStateDistance] = {Method -> "Fidelity"};

QuantumStateDistance[qs1_QuantumFiniteDimensionalState, qs2_QuantumFiniteDimensionalState, OptionsPattern[]] := If[
	!sameSizeStatesQ[qs1, qs2], $Failed,
	qdistHelper[qs1, qs2, OptionValue[Method]]];
	
QuantumStateDistance[expr___] := 
	(Message[QuantumStateDistance::inputs];$Failed);
	
Clear[sameSizeStatesQ];

sameSizeStatesQ[qs1_, qs2_] := Which[
	qs1["QuditDimension"] != qs2["QuditDimension"],
	(Message[QuantumStateDistance::dims,qs1,qs2];False),
	qs1["NumberOfQudits"] != qs2["NumberOfQudits"],
	(Message[QuantumStateDistance::numObjs,qs1,qs2];False),
	True,
	True];
	
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


(* ::Section:: *)
(*QuantumCircuit*)


Clear[QuantumCircuit];

(* Extract Properties of the Circuit *)
qcirc_QuantumCircuit[prop_String?StringQ] := extractQcircProperty[qcirc, prop];

qcirc:QuantumCircuit["BooleanFunction" -> bf_, OptionsPattern[]] := 
	If[booleanFunctionQ[bf], 
		iQCircBool[bf],
		HoldForm[qcirc]];

qcirc:QuantumCircuit[qops_List, OptionsPattern[]] :=  iQCircOps[qops];


QuantumCircuit/: MakeBoxes[qcirc1:QuantumCircuit[specs_Association, OptionsPattern[]], fmt_]:=
	circVisualize[qcirc1, fmt];


(* Converting boolean function to quantum circuit *)
boolFuncs = {Xor, And, Not, Or, Nor, Nand, Xnor};

Clear[booleanFunctionQ, booleanExpressionQ];
booleanFunctionQ[expr_] := MemberQ[boolFuncs, Head[expr]]; 
booleanExpressionQ[expr_?booleanFunctionQ] := Module[
	{subExprs},
	subExprs = expr;
	subExprs[[0]] = List;
	AllTrue[subExprs, Head[#] === Symbol || booleanFunctionQ[#]&]
	];
booleanExpressionQ[expr___] := False;


Clear[xorQ, andQ, notQ];
xorQ[expr_] := If[Head[expr] === Xor, True, False];
andQ[expr_] := If[Head[expr] === And, True, False];
notQ[expr_] := If[Head[expr] === Not, True, False];

Clear[baseLengthQ];
baseLengthQ[expr_?xorQ] := If[Length[expr] === 2, True, False];
baseLengthQ[expr_?andQ] := If[Length[expr] === 2, True, False];
baseLengthQ[expr_?notQ] := If[Length[expr] === 1, True, False];

Clear[bfReduce];

bfReduce[Not[expr_]] := {"Not", bfReduce[expr]};
bfReduce[Xor[x_, y_]] := {"Xor", bfReduce[x], bfReduce[y]};
bfReduce[expr_?xorQ] := With[{x = expr[[1]]},
	{"Xor", bfReduce[x], bfReduce[Delete[expr, 1]]}];
bfReduce[expr_?andQ] := With[{x = expr[[1]]},
	{"And", bfReduce[x], bfReduce[Delete[expr, 1]]}];
bfReduce[expr_Symbol] := expr;

Clear[esopFormQ];
esopFormQ[bf_] := With[{tf = TreeForm[bf]},
	If[
		Count[tf, _Nand, Infinity] === 0 &&
		Count[tf, _Or, Infinity] === 0 &&
		Count[tf, _Nor, Infinity] === 0 &&
		Count[tf, _Xnand, Infinity] === 0,
		True, False]
	];

Clear[convertToNormalForm];
convertToNormalForm[bf_] := With[{
	simpBF = Simplify @ BooleanConvert[bf, "ESOP"]},
	If[esopFormQ[simpBF], 
		bfReduce[simpBF],
		bfReduce @ BooleanConvert[Simplify @ bf, "ESOP"]]];

ClearAll[countANDs, countXORs, numAncillas];
countANDs[bf_] := Count[bf, "And", Infinity];
countXORs[bf_] := Count[bf, "Xor", Infinity];
numAncillas[boolFun_] := countANDs[boolFun] + countXORs[boolFun] -1;


Clear[extractQcircProperty];
extractQcircProperty[qcirc_, "QuditDimension"] :=  qcirc[[1, "QuditDimension"]];
extractQcircProperty[qcirc_, "Instructions"] :=  qcirc[[1, "Instructions"]];
extractQcircProperty[qcirc_, "Arity"] :=  qcirc[[1, "Arity"]];
extractQcircProperty[qcirc_, "NumberOfAncillas"] :=  qcirc[[1, "NumberOfAncillas"]];
extractQcircProperty[qcirc_, "CircuitDiagram"] := circuitRepresentation[qcirc];


Clear[iQCircOps];
iQCircOps[ops_] := Module[
	{wires, arity, dim, assoc},
	wires = Range @@ MinMax[Union[ops[[All, 2]]]];
	arity = Length[wires];
	dim = ops[[1,1]][[1,"QuditDimension"]];
	assoc = <|{
		"Instructions" -> ops,
		"QuditDimension" -> dim,
		"Arity" -> arity,
		"Wires" -> wires,
		"NumberOfAncillas" -> 0,
		"InitializeAncillas" -> False,
		"TraceOutAncillas" -> False
		}|>;
	QuantumCircuit[assoc]
	];

				
iQCircBool[bf_] := Module[
	{vars, normalBF, nAnc, ancCount, instructions, ANDPos, XORPos, assoc},
	vars = BooleanVariables[bf];
	normalBF = convertToNormalForm[bf];
	nAnc = numAncillas[normalBF];
	ANDPos = Map[Drop[#, -1]&, Reverse[SortBy[Position[normalBF, "And"], Depth]]];
	XORPos = Map[Drop[#, -1]&, Reverse[SortBy[Position[normalBF, "Xor"], Depth]]];
	
	{ancCount, normalBF, instructions} = 
		Fold[boolCircANDHelper[#1, #2, nAnc]&, {0, normalBF, {}}, ANDPos];
	{ancCount, normalBF, instructions} = 
		Fold[boolCircXORHelper[#1, #2, nAnc]&, {ancCount, normalBF, instructions}, XORPos];	
	instructions = boolCircNegationHelper[vars, instructions];
	instructions = cleanGarbage[instructions];
	(* Turn instructions into matrix operations *)
	instructions = instructionsToOps[instructions, Sort @ vars, nAnc];
	assoc = <|{
		"Instructions" -> instructions,
		"QuditDimension" -> 2,
		"Arity" -> Length[vars],
		"NumberOfAncillas" -> nAnc,
		"InitializeAncillas" -> True,
		"TraceOutAncillas" -> True,
		"Measure" -> False,
		"QuantumResultRegister" -> True,
		"ClassicalResultRegister" -> False,
		"Wires" -> Join[Sort[vars], {"resReg"}, 
			Map[StringJoin["anc[", ToString[#], "]"]&, Range[1, nAnc]]]
		}|>;
	
	QuantumCircuit[assoc]
	];
	
ClearAll[instructionToOp, instructionsToOps, instrToOpHelper];
instructionsToOps[instrs_, vars_, nAnc_] := Map[instructionToOp[#, vars, nAnc]&, instrs];

instructionToOp[instr_, vars_, nAnc_] := Which[
	instr[[1]] === "Toffoli", 
	QuantumMatrixOperation["Toffoli"] -> 
		{instrToOpHelper[instr[[2]], vars, nAnc],
		 instrToOpHelper[instr[[3]], vars, nAnc],
		 instrToOpHelper[instr[[4]], vars, nAnc]},
	instr[[1]] === "CNOT",
	QuantumMatrixOperation["CNOT"] -> 
		{instrToOpHelper[instr[[2]], vars, nAnc],
		 instrToOpHelper[instr[[3]], vars, nAnc]},
	instr[[1]] === "SigmaX",
	QuantumMatrixOperation["SigmaX"] -> 
		{instrToOpHelper[instr[[2]], vars, nAnc]}
		];

instrToOpHelper[var_, vars_, nAnc_] := Which[
	MemberQ[vars, var],
	Position[vars, var][[1,1]],
	var === "ResultRegister",
	Length[vars] + 1,
	var[[1]] === "Ancilla",
	var[[2]] + Length[vars] + 2];


Clear[boolCircANDHelper];
boolCircANDHelper[{anc_, bf_, instr_}, ANDpos_, nAnc_] := Module[
	{and, args},
	and = Extract[bf, ANDpos];
	args = If[ANDpos === {}, 
		{bf[[2]], bf[[3]]},
		{and[[2]], and[[3]]}];
	{anc + 1, ReplacePart[bf, ANDpos -> ("Ancilla" -> anc)], 
		Join[instr, {{"Toffoli", args[[1]], args[[2]], 
			If[nAnc === anc, "ResultRegister", ("Ancilla" -> anc)]}}]} 
	];
	
Clear[boolCircXORHelper];
boolCircXORHelper[{anc_, bf_, instr_}, XORpos_, nAnc_] := Module[
	{xor, args},
	xor = Extract[bf, XORpos];
	args = If[XORpos === {}, 
		{bf[[2]], bf[[3]]},
		{xor[[2]], xor[[3]]}];
	{anc + 1, ReplacePart[bf, XORpos -> ("Ancilla" -> anc)], 
		Join[instr, {{"CNOT", args[[1]], If[nAnc === anc, "ResultRegister", ("Ancilla" -> anc)]},
			{"CNOT", args[[2]], If[nAnc === anc, "ResultRegister", ("Ancilla" -> anc)]}}]}
	];
	
ClearAll[boolCircNegationHelper, negationHelper];
boolCircNegationHelper[vars_, instr_] := Module[
	{negVars, notPos, newInstr},
	negVars = Association @ Map[# -> False &, vars];
	notPos = Reverse @ Position[instr, "Not"];
	{negVars, newInstr} = Fold[negationHelper[#1, #2]&, {negVars, instr}, notPos];
	newInstr
	];
	
negationHelper[{negVars_, instr_}, NOTpos_] := Module[
	{var, newInstr, newNegVars},
	If[instr[[NOTpos[[1]], 1]] === "SigmaX",
		Return[negationHelper[{negVars, instr}, {NOTpos[[1]] + 1, NOTpos[[2]], NOTpos[[3]]}]]];
	var = Extract[instr, {NOTpos[[1]], NOTpos[[2]], 2}];
	newInstr = ReplacePart[instr, {NOTpos[[1]], NOTpos[[2]]} -> var];
	newNegVars = negVars;
	If[negVars[var] === False,
		newInstr = Insert[newInstr, {"SigmaX", var}, NOTpos[[1]]]];
	newNegVars[var] = True;
	{newNegVars, newInstr}
	];

Clear[cleanGarbage];
cleanGarbage[instr_] := With[{ndrop = Count[instr, "ResultRegister", Infinity]},
	Join[instr, Reverse[Drop[instr, - ndrop]]]];


(* Front End Fornatting of Circuit *)
Clear[circVisualize];
Clear[circVisualize];
circVisualize[circ_QuantumCircuit, fmt_] := 
	BoxForm`ArrangeSummaryBox[QuantumCircuit,
		circLogo,
		circLogo,
		circBaseVisual[circ], 
		circExpandedVisual[circ], 
		fmt];
		
Clear[circBaseVisual, circExpandedVisual];
ClearAll[circArityVisualize, circNumAncillasVisualize];

circBaseVisual[circ_] := 
	Join[circArityVisualize[circ], circNumAncillasVisualize[circ]];

circExpandedVisual[circ_] := Join[
	circDimVisualize[circ], {}];
	
circDimVisualize[circ_] := With[{
	dim = circ[[1,"QuditDimension"]]},
	If[! dim ===  2,
		{BoxForm`MakeSummaryItem[
			{"Qudit Dimension: ", dim},
			 StandardForm]},
		{}
		]
	]
	
circArityVisualize[circ_] := 
	{BoxForm`MakeSummaryItem[
		{"Arity: ",circ[[1,"Arity"]]},
		 StandardForm]};
		 
circNumAncillasVisualize[circ_] :=
	{BoxForm`MakeSummaryItem[
		{"Number of Ancillas: ", circ[[1,"NumberOfAncillas"]]},
		 StandardForm]};





baseString[op_] := If[MemberQ[Keys[op[[1]]], "OperationName"],
	op[[1, "OperationName"]],
	"\[ScriptCapitalU]"];



(*  Circuit Representation Graphics Generation *)
wireGraphics[numSteps_, wires_] := Flatten[MapIndexed[
		{{Text[ToString[#1], {-3 * numSteps, 100 * First[#2]}]},
		{AbsoluteThickness[2], 
		Line[{{10,100* First[#2]}, {60*numSteps + 60,100 *First[#2]}}]}}&, 
		Reverse @ wires]];
		
Clear[circLogo];
circLogo := Graphics[Join[
			{EdgeForm[GrayLevel[0]],FaceForm[], Rectangle[{-30,0},{870, 400}]},
				wireGraphics[12, {"", "", ""}]], ImageSize->{50,30}];
				
				
Clear[initialWirePositions];
initialWirePositions[wires_] := Association @ Thread[wires -> 0];

Clear[updateWirePositions];
updateWirePositions[posAssoc_, w_, newPos_] := Module[
	{newAssoc},
	newAssoc = posAssoc;
	newAssoc[w] = newPos;
	newAssoc
];

Clear[addOpToCircGraphicsSteps];
addOpToCircGraphicsSteps[{numSteps_, positions_, steps_}, op_] := Module[
	{wsRange, pos, ws, rightSteps, offset, newPositions, newSteps},
	ws = op[[2]];
	wsRange = Range @@ MinMax[ws];
	pos = Max[Map[positions, wsRange]];
	newPositions = Fold[updateWirePositions[#1, #2, pos + 1]&, positions, wsRange];
	If[pos === numSteps, 
		newSteps = Append[steps, {op}]; Return[{numSteps + 1, newPositions, newSteps}]];
	newSteps = steps;
	newSteps[[pos + 1]] = Join[steps[[pos + 1]], {op}];
	{numSteps, newPositions, newSteps}
];	

Clear[circuitGraphicsSteps];
circuitGraphicsSteps[circ_] := Module[
	{ops, positions, numSteps, steps, wires},
	ops = circ[[1,"Instructions"]];
	wires = Range @ Length[circ[[1,"Wires"]]];
	ops = Map[#[[1]] -> (Length[wires] + 1 - #[[2]])&, ops];
	positions = initialWirePositions[wires];
	numSteps = 0;
	steps = {};
	{numSteps, positions, steps} = 
		Fold[addOpToCircGraphicsSteps[#1, #2]&, {numSteps, positions, steps}, ops];
	steps
];

	
Clear[cnotGraphics];
cnotGraphics[wireNums_, stepNum_] := Module[
	{cntrl, trgt, rad, lineSeg},
	cntrl = wireNums[[1]];
	trgt = wireNums[[2]];
	rad = 20;
	lineSeg = If[cntrl > trgt, - rad, rad];
	{
	PointSize[.07],
	Point[{60*stepNum, 100*cntrl}], Black,
	Circle[{60*stepNum, 100*trgt}, rad],
	Thick,Black,Line[{{60*stepNum, 100*trgt + lineSeg},{60*stepNum, 100*cntrl}}]
	}
	];
	
Clear[cphaseGraphics];
cphaseGraphics[wireNums_, stepNum_] := Module[
	{cntrl, trgt, connector, box},
	cntrl = wireNums[[1]];
	trgt = wireNums[[2]];
	connector = Line[{{60*stepNum, 100*trgt},{60*stepNum, 100*cntrl}}];
	box = {
		EdgeForm[Black],
		FaceForm[LightBlue], 
		Rectangle[
			{60 * stepNum -15, 100 * trgt -20},
			{60 * stepNum + 15, 100 * trgt + 20}
			],
		Text["\[ScriptCapitalZ]",{60*stepNum, 100*trgt}]
		};
	{
		PointSize[.07], 
		Point[{60*stepNum, 100*cntrl}],
		Thick, Black, 
		connector,
		box
	}
	];
	
Clear[toffGraphics];
toffGraphics[wireNums_, stepNum_] := Module[
	{cntrl1, cntrl2, trgt, rad, lowLine, highLine, min, max},
	cntrl1 = wireNums[[1]];
	cntrl2 = wireNums[[2]];
	trgt = wireNums[[3]];
	{min, max} = MinMax[wireNums];
	
	rad = 20;
	lowLine = If[trgt === min, -rad, 0];
	highLine = If[trgt === max, rad, 0];
	{
	PointSize[.07],
	Point[{60*stepNum, 100*cntrl1}], 
	Point[{60*stepNum, 100*cntrl2}], Black,
	Circle[{60*stepNum, 100*trgt}, rad],
	Thick,Black,Line[{{60*stepNum, 100*min + lowLine},{60*stepNum, 100*max + highLine}}]
	}
	];
	
Clear[fredkinGraphics];
fredkinGraphics[wireNums_, stepNum_] := Module[
	{cntrl, sw1, sw2, min, max, l1, l2, l3, l4, connector},
	cntrl = wireNums[[1]];
	sw1 = wireNums[[2]];
	sw2 = wireNums[[3]];
	{min, max} = MinMax[wireNums];
	
	l1 = Line[{{60*stepNum - 10, 100*sw1 - 15},{60*stepNum + 10, 100*sw1 + 15}}];
	l2 = Line[{{60*stepNum - 10, 100*sw1 + 15},{60*stepNum + 10, 100*sw1 - 15}}];
	l3 = Line[{{60*stepNum - 10, 100*sw2 - 15},{60*stepNum + 10, 100*sw2 + 15}}];
	l4 = Line[{{60*stepNum - 10, 100*sw2 + 15},{60*stepNum + 10, 100*sw2 - 15}}];
	connector = Line[{{60*stepNum, 100*min},{60*stepNum, 100*max}}];
	{
		PointSize[.07], 
		Point[{60*stepNum, 100*cntrl}], 
		Thick, Black, 
		l1, l2, l3, l4, 
		connector
	}
	];
	
Clear[deutschGraphics];
deutschGraphics[wireNums_, stepNum_] := Module[
	{cntrl1, cntrl2, trgt, min, max, box, connector, rot},
	cntrl1 = wireNums[[1]];
	cntrl2 = wireNums[[2]];
	trgt = wireNums[[3]];
	{min, max} = MinMax[wireNums];
	
	box = {
		EdgeForm[Black],
		FaceForm[LightBlue], 
		Rectangle[
			{60 * stepNum -25, 100 * trgt -20},
			{60 * stepNum + 25, 100 * trgt + 20}
			],
		Text["R(\[Theta])",{60*stepNum, 100*trgt}]
		};
	
	connector = Line[{{60*stepNum, 100*min},{60*stepNum, 100*max}}];
	{
		PointSize[.07], 
		Point[{60*stepNum, 100*cntrl1}],
		Point[{60*stepNum, 100*cntrl2}], 
		Thick, Black, 
		connector,
		box
	}
	];
	
Clear[swapGraphics];
swapGraphics[wireNums_, stepNum_] := Module[
	{l1, l2, l3, l4, connector, w1, w2},
	w1 = wireNums[[1]];
	w2 = wireNums[[2]];
	l1 = Line[{{60*stepNum - 10, 100*w1 - 15},{60*stepNum + 10, 100*w1 + 15}}];
	l2 = Line[{{60*stepNum - 10, 100*w1 + 15},{60*stepNum + 10, 100*w1 - 15}}];
	l3 = Line[{{60*stepNum - 10, 100*w2 - 15},{60*stepNum + 10, 100*w2 + 15}}];
	l4 = Line[{{60*stepNum - 10, 100*w2 + 15},{60*stepNum + 10, 100*w2 - 15}}];
	connector = Line[{{60*stepNum, 100*w1},{60*stepNum, 100*w2}}];
	{Thick, Black, l1, l2, l3, l4, connector}
	];
	
Clear[rootSwapGraphics];
rootSwapGraphics[wireNums_, stepNum_] := Module[
	{l1, l2, l3, l4, connector1, connector2, w1, w2, circle, min, max, txt},
	w1 = wireNums[[1]];
	w2 = wireNums[[2]];
	{min, max} = MinMax[{w1, w2}];
	l1 = Line[{{60*stepNum - 10, 100*w1 - 15},{60*stepNum + 10, 100*w1 + 15}}];
	l2 = Line[{{60*stepNum - 10, 100*w1 + 15},{60*stepNum + 10, 100*w1 - 15}}];
	l3 = Line[{{60*stepNum - 10, 100*w2 - 15},{60*stepNum + 10, 100*w2 + 15}}];
	l4 = Line[{{60*stepNum - 10, 100*w2 + 15},{60*stepNum + 10, 100*w2 - 15}}];
	connector1 = Line[{{60*stepNum, 100*min},{60*stepNum, 50*(min+max)-20}}];
	connector2 = Line[{{60*stepNum, 100*max},{60*stepNum, 50*(min+max)+20}}];
	circle = Circle[{60*stepNum, 50*(w1+w2)}, 20];
	txt = Text["1/2",{60*stepNum, 50*(w1+w2)}];
	{Thick, Black, l1, l2, l3, l4, connector1, connector2, circle, txt}
	];
	
Clear[circuitRepresentation];
circuitRepresentation[circ_] := Module[
	{test, wires, steps, width, height, imageW, imageH, xScroll, yScroll, wireLines, graphics},
	wires = circ[[1,"Wires"]];
	steps = circuitGraphicsSteps[circ];
	width = circuitWidth[steps];
	height = circuitHeight[wires];
	imageW = circuitImageWidth[steps];
	imageH = circuitImageHeight[wires];
	
	xScroll = circuitXScroll[steps];
	yScroll = circuitYScroll[wires];
	
	wireLines = wireGraphics[Length[steps], wires];
	graphics = Flatten[ 
		MapIndexed[{#1, wires, #2[[1]]}&, steps, {2}], 1];
	graphics = Map[opGraphics, graphics];
	Pane[Panel @ Panel @ Graphics[Join[wireLines, graphics], ImageSize->{imageW, imageH}], ImageSize->{400,200}, Scrollbars->{xScroll, yScroll}]
	];
	
ClearAll[circuitWidth, circuitHeight];
circuitWidth[steps1_] := 60 + 60 * Length[steps1];
circuitHeight[wires_] := 100 * Length[wires];

ClearAll[circuitImageWidth, circuitImageHeight];
circuitImageWidth[stps_] := Min[300, circuitWidth[stps]];
circuitImageHeight[wires_] := Min[100, circuitHeight[wires]];

ClearAll[circuitXScroll, circuitYScroll];
circuitXScroll[steps2_] := If[Length[steps2] >= 5, True, False];
circuitYScroll[wires_] := If[Length[wires] >= 5, True, False];

Clear[opGraphics];
opGraphics[{operation_ -> wires_, allWires_, stepNum_}] := 
	If[Length[wires] > 1, 
		multiQuditOpGraphics[wires, operation, stepNum, Length @ allWires],
		singleQuditOpGraphics[wires[[1]], operation, stepNum, Length @ allWires]
		];
	
			
Clear[singleQuditOpGraphics];
singleQuditOpGraphics[wire_, operation_, stepNum_, numWires_] := 
	{EdgeForm[Black],FaceForm[LightBlue], Rectangle[{60 * stepNum -15, 100 * wire -20},{60 * stepNum + 15, 100 * wire + 20}],
	Text[singleQuditGraphicsHelper[baseString[operation]],{60*stepNum, 100*wire}]};	

	
Clear[singleQuditGraphicsHelper];
singleQuditGraphicsHelper["SigmaX"] := "\!\(\*SubscriptBox[\(\[Sigma]\), \(\[ScriptCapitalX]\)]\)";
singleQuditGraphicsHelper["SigmaY"] := "\!\(\*SubscriptBox[\(\[Sigma]\), \(\[ScriptCapitalY]\)]\)";
singleQuditGraphicsHelper["SigmaZ"] := "\!\(\*SubscriptBox[\(\[Sigma]\), \(\[ScriptCapitalZ]\)]\)";
singleQuditGraphicsHelper["Hadamard"] := "\[ScriptCapitalH]";
singleQuditGraphicsHelper[expr___] := expr;
	
			
Clear[multiQuditOpGraphics];
multiQuditOpGraphics[wires_, operation_, stepNum_, numWires_] := 
	With[
	{base = baseString[operation]},
	(Which[
		base === "CNOT",
		cnotGraphics[wires, stepNum],
		
		base === "SWAP",
		swapGraphics[wires, stepNum],
		
		base === "RootSWAP",
		rootSwapGraphics[wires, stepNum],
		
		base === "CPHASE",
		cphaseGraphics[wires, stepNum],
		
		base === "Toffoli",
		toffGraphics[wires, stepNum],
		
		base === "Fredkin",
		fredkinGraphics[wires, stepNum],
		
		ListQ[base] && base[[1]] === "Deutsch",
		deutschGraphics[wires, stepNum],
		
		True,
		"ERROR"
		])
	];


(* ::Section:: *)
(*QuantumEvaluate*)


Clear[QuantumEvaluate];

(* Options *)
Options[QuantumEvaluate] = {"Trials" -> 1, "BackEnd" -> "ClassicalSimulator"};

QuantumEvaluate /: qev:QuantumEvaluate[order_ -> qop_QuantumMatrixOperation, qs_QuantumFiniteDimensionalState, OptionsPattern[]]:=
	If[
		explicitQuantumMatrixOpQ[qop] &&
		explicitQuantumStateQ[qs],
		iQOpEval[qop, qs, order, OptionValue["BackEnd"]],
		HoldForm[qev]
		];
		
QuantumEvaluate /: qev:QuantumEvaluate[order_ -> qmeas_QuantumMeasurement, qs_QuantumFiniteDimensionalState, OptionsPattern[]]:=
	If[
		explicitQuantumMeasurementQ[qmeas] &&
		explicitQuantumStateQ[qs],
		iQMeasEval[qmeas, qs, order, OptionValue["Trials"], OptionValue["BackEnd"]],
		HoldForm[qev]
		];
		
Clear[iQOpEval];
iQOpEval[qop_, qs_, order_, backend_] := 
	If[backend === "ClassicalSimulator",
		iQOpSimulate[qop, qs, order],
		iQOpExecute[qop, qs, order]];
		
Clear[iQMeasEval];
iQMeasEval[qmeas_, qs_, order_, trials_, backend_] := 
	If[backend === "ClassicalSimulator",
		iQMeasSimulate[qmeas, qs, order, trials],
		iQMeasExecute[qmeas, qs, order, trials]];
		
Clear[iQCircEval];
iQCircEval[qcirc_, qs_, order_, trials_, backend_] := 
	If[backend === "ClassicalSimulator",
		iQCircSimulate[qcirc, qs, order, trials],
		iQCircExecute[qcirc, qs, order, trials]];
		


Clear[iQOpSimulate];
iQOpSimulate[qop_, qs_, order_] := With[{
	mat = If[qop["Arity"] === 1, 
		singleQuditOpMat[qop, qs, order[[1]]],
		multiQuditOpMat[qop, qs, order]]},
		applyMatToState[mat, qs]];
	
Clear[iQMeasSimulate];
iQMeasSimulate[qmeas_, qs_, order_, trials_] := If[
	qmeas["MeasurementType"] === "Projection",
	iQProjectionSimulate[qmeas, qs, order, trials],
	iQPOVMSimulate[qmeas, qs, order, trials]];	
	
Clear[eigToProjector];
eigToProjector[eig_] := SparseArray @ ConjugateTranspose[{eig}].{eig}
		
Clear[iQProjectionSimulate];
iQProjectionSimulate[qmeas_, qs_, order_, trials_] := With[{
	projectors = If[qmeas["Arity"] === 1,
	singleQuditProjector[qmeas, qs, order[[1]]],
	multiQuditProjector[qmeas, qs, order]]},
	applyProjectorsToState[projectors, qs, qmeas, trials]];

	
Clear[applyProjectorsToState];
applyProjectorsToState[projectors_, qs_, qmeas_, trials_] := Module[
	{dm, probs},
	dm = qs["DensityMatrix"];
	probs = Map[Tr[#.dm]&, projectors];
	RandomChoice[probs -> qmeas["Eigenvalues"], trials]
	];


Clear[iQPOVMSimulate];
iQPOVMSimulate[qmeas_, qs_, order_, trials_] := With[{
	povm = If[qmeas["Arity"] === 1,
	singleQuditPOVM[qmeas, qs, order[[1]]],
	multiQuditPOVM[qmeas, qs, order]]},
	applyPOVMToState[povm, qs, qmeas, trials]];
	
	
Clear[applyPOVMToState];
applyPOVMToState[povm_, qs_, qmeas_, trials_] := Module[
	{dm, probs},
	dm = qs["DensityMatrix"];
	probs = Map[Tr[#.dm]&, povm];
	RandomChoice[probs -> Range[Length[probs]], trials]
	];


						
Clear[applyMatToState];
applyMatToState[mat_, state_] := Module[
	{newState},
	newState = state;
	If[pureStateQ[state], 
		newState[[1,"StateVector"]] = FullSimplify @ SparseArray[mat.state["StateVector"]], 
		newState[[1,"DensityMatrix"]] = FullSimplify @  SparseArray[mat.state["DensityMatrix"].ConjugateTranspose[mat]]
		];
	newState
	];
	
	
	
	
	
ClearAll[singleQuditOpMat, multiQuditOpMat];
singleQuditOpMat[qop_, qs_, pos_] := Module[
	{opMats, dim},
	dim = qs["QuditDimension"];
	opMats = ConstantArray[id[dim], qs["NumberOfQudits"]];
	opMats[[pos]] = qop["MatrixRepresentation"];
	If[Length @ opMats == 1, opMats[[1]], KroneckerProduct @@ opMats]
	];
	
multiQuditOpMat[qop_, qs_, pos_] := Module[
	{qObjs, arity, n, dim, passiveObjs, order, passiveRep, shape, opMat, tpLevels},
	n = qs["NumberOfQudits"];
	qObjs = Range[n];
	arity = qop["Arity"];
	dim = qs["QuditDimension"];
	passiveObjs = Complement[qObjs, pos];
	order = Join[pos, passiveObjs];
	passiveRep = ConstantArray[id[dim], Length @ passiveObjs];
	shape = ConstantArray[dim, 2 * n];
	opMat = If[Length @ passiveObjs === 0,
		qop["MatrixRepresentation"],
		KroneckerProduct @@ Join[{qop["MatrixRepresentation"]}, passiveRep]];
	opMat = ArrayReshape[opMat, shape];
	tpLevels = InversePermutation[Ordering @ order];
	tpLevels = Join[tpLevels, tpLevels + n];
	opMat = Transpose[opMat, tpLevels];
	opMat = ArrayReshape[Flatten @ opMat, {dim^n, dim^n}];
	opMat
	];


End[];
EndPackage[];
