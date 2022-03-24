# NMRSpectraSimulator

[![Build Status](https://github.com/RoyCCWang/NMRSpectraSimulator.jl/workflows/CI/badge.svg)](https://github.com/RoyCCWang/NMRSpectraSimulator.jl/actions)
[![Coverage](https://codecov.io/gh/RoyCCWang/NMRSpectraSimulator.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/RoyCCWang/NMRSpectraSimulator.jl)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ai4dbiological-systems.github.io/NMRSpectraSimulator.jl/)
Uses the spin Hamiltonian theory for isotropic liquid 1D 1H NMR spectroscopy. See Ernst's "Principles of Nuclear Magnetic Resonance in One and Two Dimensions" (ISBN 0-19-855647-0) and Levitt's "Spin Dynamics" (ISBN 978-0-470-51117-6) for the theory.


## To Install:
Since the pre-requisite package GISSMOReader.jl is an unregistered with the Julia public registry, you must manually install them first. Open a Julia REPL session, and run the following commands.
```
import Pkg
Pkg.add(path="https://github.com/AI4DBiological-Systems/GISSMOReader.jl")
Pkg.add(path="https://github.com/AI4DBiological-Systems/NMRSpectraSimulator.jl")
```

 
## Examples (in progress)
The preliminary library description and tutorial (in progress) website is located [here](https://ai4dbiological-systems.github.io/NMRSpectraSimulator.jl/)

More examples in the example folder. We will continue to use Weave.jl to work off the existing example scripts to finish annoted tutorials.
