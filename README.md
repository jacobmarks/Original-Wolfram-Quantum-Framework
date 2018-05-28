# QuantumComputing

Symbolic Quantum Computing Package for the Wolfram Language. Includes support for qudits, mixed states, and generation of reversible quantum circuits from classical boolean functions. This package began as a project during the Wolfram Summer School 2017, and grew into an internship in Algorithms R&D for Wolfram Research.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Installing

Download a copy of QuantumComputing.m and place this in the directory with the notebook you are trying to run. In an evaluation cell in your notebook, paste and evaluate the following command:

```
Get[FileNameJoin[NotebookDirectory[], "QuantumComputing.m"]]
```


## Capabilities

The symbol notebooks illustrate the full scope of capabilities of the package. The documentation includes usage syntax, examples, and more involved applications. This section does not attempt to duplicate the complete documentation for the package. Instead, its purpose is to introduce the reader to the essentials.

### States

`QuantumFiniteDimensionalState` is the generic symbol used for declaring quantum states. Continuous states are not yet supported. Create a two-level quantum system, or "qubit", in the "0" computational basis state:

```
qs = QuantumFiniteDimensionalState[{{"BasisState", {0}}]
```
Extract information about the quantum state such as its purity, Von Neumann entropy, and number of qudits using keywords. For example, 

```
qs["BasisStates"]
```

returns a list of the canonical basis states for the system. Keywords include

* StateVector
* DensityMatrix
* Purity
* VonNeumannEntropy
* BasisStates
* Plot
* BlochPlot
* NumberOfQudits
* QuditDimension
* PureStateQ
* MixedStateQ

QuantumFiniteDimensionalState accepts as inputs both numerical (or symbolic) arrays and keywords for recognized states. To make a symbolic qutrit with complex coefficients a, b and c for basis states 0, 1, and 2 respectively, one would write

```
QuantumFiniteDimensionalState[{{a, b, c}}, "QuditDimension" -> 3]
```

Note qudits of arbitrary dimension, or spin-d/2 particles, are fully supported. This support extends to the measurements and matrix operations described below.

Alternatively, one could generate the PhiPlus Bell state with

```
QuantumFiniteDimensionalState[{"PhiPlus"}]
```

Some recognized states include

* Plus
* Minus
* Left
* Right
* PhiPlus
* PhiMinus
* PsiPlus
* PsiMinus
* GHZ
* W

New quantum states can be obtained by modifying existing quantum states using the symbols `QuantumProduct`,`QuantumProduct`, and `QuantumPartialTr`. `QuantumProduct` generates tensor product states. `QuantumMixture` generates a mixed quantum state from a statistical ensemble of quantum states with given classical probabilities. And `QuantumPartialTr` can be used to trace out speficied subsystems from a quantum state. 

### Operators

Operators are split into two classes: `QuantumMatrixOperation`, which includes but is not limited to the unitary operations that serve as gates in a quantum circuit, and `QuantumMeasurement`, which can be used to specify either a projective measurement, traditionally representing an observable, or a generalized positive operator valued measurement (POVM).

For matrix operations, for example, one could generate a Hadamard gate either explicitly via input array:

```
h = 1/Sqrt[2] {{1, 1}, {1, -1}};
QuantumMatrixOperation[h]
```
or via keyword:

```
QuantumMatrixOperation["Hadamard"]
```

With the "QuditDimension" Option, one can obtain the d-level version of a quantum matrix operation if it is sensible. For instance:

```
QuantumMatrixOperation["SigmaZ", "QuditDimension" -> 5]
```

generates the Pauli Z operator acting on 5-level systems.

Recognized matrix operations include

* SigmaX
* SigmaY
* SigmaZ
* SigmaPlus
* SigmaMinus
* Hadamard
* S
* T
* CNOT
* CPHASE
* SWAP
* Deutsch
* Fredkin
* Toffoli
* Fourier
* SUM
* RandomUnitary
* RotX
* RotY
* RotZ

For measurements, you must specify the type of measurement as either "Observable" or "POVM". For example, to instantiate a projective measurement in the x basis, one would call

```
QuantumMatrixOperation["Observable" -> "SigmaX"]
```

Like quantum states, both matrix operations and measurements have properties that can be extracted.

For `QuantumMatrixOperation`, these include

* MatrixRepresentation
* Arity
* QuditDimension
* UnitaryQ
* HermitianQ

and for `QuantumMeasurement`,

* Arity
* QuditDimension
* MeasurementType
* POVM (POVM only)
* Eigenvalues (Observable only)
* Eigenvectors (Observable only)

### Circuits

In `QuantumComputing`, circuits are represented by the `QuantumCircuit` symbol. These circuits consist of unitary operations, projective measurements, and control flow that allows for classically conditioned operations and repeat-until-success (RUS) measurements, which allow for the simulation of nonlinearity with finite dimensional quantum states. 

circuits can be concatenated effortlessly with "+", and reversible quantum circuits for (quantum versions of) Boolean functions can be generated automatically with

```
QuantumCircuit["BooleanFunction" -> bf]
```

where bf is a BooleanFunction in Mathematica.

## Authors

* **Jacob Marks**  - (jamarks@stanford.edu or jacobm@wolfram.com)


## Acknowledgments

* Stephen Wolfram
* Jose Martin-Garcia
* Cesar Guerra
* Xavier Roy
* Wolfram Research Algorithms R&D
