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

### QuantumFiniteDimensionalState

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
* PureStateQ
* MixedStateQ


## Authors

* **Jacob Marks**  - (jamarks@stanford.edu or jacobm@wolfram.com)


## Acknowledgments

* Stephen Wolfram
* Jose Martin-Garcia
* Cesar Guerra
* Xavier Roy
