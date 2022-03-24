"""
title: NMRSpectraSimulator.jl Documentation
author Roy Wang
date: 2022-03-23
"""

# Library description

Given a list of user-specified compounds and parameters, this library provides the simulated resonance components of the 1D 1H NMR experiment signal. The signal model is taken to be the complex-valued free-induction decay (FID) model 

$f_{\text{FID}}\left(t,\,w,\,\alpha,\,\Omega,\,\lambda,\,\beta\right)=\sum_{n,l_n}w_{n}\alpha_{n,l_n}e^{i\left(\Omega_{n,l_n}t-\beta_{n,l_n}\right)}e^{-\lambda_{n,l_n}t},\,t>0,$
and its Fourier transform, called the FID spectra model, is 

Here, $n$ is the compound index, $l_n$ is the resonance component index for compound $n$. The parameters of this signal model are: The compound molar concentration $w$, the resonance frequency $\Omega$, the resonance amplitude $\alpha$, the resonance T2 decay $\lambda$, the resonance phase $\beta$. The resonance frequencies and amplitudes of a compound can be simulated if the physical chemistry parameters, chemical shift and J-coupling, are known. This library uses the physical chemistry parameters from the [GISSMO website](http://gissmo.nmrfam.wisc.edu/library), which were obtained by empirically fitting a different forward model than the one used in this library against real NMR experiments. For GISSMO-related details, see [Dashti et. al 2017](https://pubs.acs.org/doi/abs/10.1021/acs.analchem.8b02660) and [Dashti et. al 2018](http://pubs.acs.org/doi/abs/10.1021/acs.analchem.7b02884) for their official journal publications.


We do not recommend allowing each resonance component to have its only degree of freedom when fitting this simulated signal to data. This is because there are a large number of resonance coherences (e.g., up to 3003 for alpha-D-Glucose). To address this application challenge, this library also provides grouping information for the resonance components, such that the changes in frequency that one sees across experiments for the components in a frequency group are approximately constant within each group. A separate grouping for changes in T2 decay and phase is also provided.

These grouping information and the resonance components are used by this library to provide fast spline-based surrogate models to the simulated spectra. This surrogate model is what is used in NMRCalibrate.jl to do model fitting against data.

# What this tutorial is about.
The workflow is 


# Getting Started

## Tutorial for simulating a mixture of compounds.

#Describe what this tutorial is about.

```julia; echo=false, results="hidden"
import NMRSpectraSimulator

using LinearAlgebra
using FFTW
import PyPlot
import Statistics

import Random
Random.seed!(25)
```

Two types of parameters: the coherence simulation and spectrometer/experiment parameters.


Coherence-related parameters for the spectrum simulation.

```julia
tol_coherence = 1e-2
α_relative_threshold = 0.05
λ0 = 3.4
Δcs_max = 0.2
κ_λ_lb = 0.5
κ_λ_ub = 2.5
```
Δcs_max


Spectrometer-related parameters. The [NMRDataSetup.jl](https://github.com/AI4DBiological-Systems/NMRDataSetup.jl) library get these parameters from user-supplied NMR experiments.
```julia
fs = 14005.602240896402 # sampling frequency, in Hz.
SW = 20.0041938620844 # spectral window, in ppm. Described in Bruker TopSpin manual,
ν_0ppm = 10656.011933076665 # where the 0 ppm reference freuency should be, in Hz.
```

NMR spectroscopy uses ppm as the frequency unit. The physics simulation use frequency units in Hz. The following are the conversion formulae:
```julia
hz2ppmfunc = uu->(uu - ν_0ppm)*SW/fs
ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)
```


We need to provide the path to the location of the [/src/input/molcules](https://github.com/AI4DBiological-Systems/NMRData/tree/main/src/input/molecules) in the [NMRData.jl](https://github.com/AI4DBiological-Systems/NMRData) repository. This folder contain the stored J-coupling and chemical shift constants from the [GISSMO library](http://gissmo.nmrfam.wisc.edu/library). Stored in JLD file format.

```julia
base_path_JLD = "/home/roy/Documents/repo/NMRData/src/input/molecules"
```
Note that the chemical shift tend to change slightly (within +/- 0.1 ppm) for every NMR experiment. The J-coupling constants for a compound shouldn't change in theory, but the values provided by the GISSMO library are empirically fitted, thus prone to misfit error and physics model descrepancy. In [NMRCalibrate.jl](https://github.com/AI4DBiological-Systems/NMRCalibrate.jl), we work on the fitting of the simulated spectra (provided in this package) to NMR experiment data.


In this tutorial, we simulate following compounds:
```julia
molecule_names = ["L-Serine"; "L-Phenylalanine"; "DSS";]
```
See [GISSMOReader.jl](https://github.com/AI4DBiological-Systems/GISSMOReader.jl) for the current list of feasible compounds.





















Provide the range of frequency shift that we 
```julia
Δcs_max_mixture = collect( Δcs_max for i = 1:length(molecule_names))


dummy_SSFID = NMRSpectraSimulator.SpinSysFIDType1(0.0)
@time mixture_params = NMRSpectraSimulator.setupmixtureproxies(molecule_names,
    base_path_JLD, Δcs_max_mixture, hz2ppmfunc, ppm2hzfunc, fs, SW, λ0,
    ν_0ppm, dummy_SSFID;
    tol_coherence = tol_coherence,
    α_relative_threshold = α_relative_threshold)
As = mixture_params
```








We fcan apply this code:

```julia
f = xx->xx^2
f(3)
```

support inline eval: $f(4) = `j f(4)`$



Do plot.

```julia; fig_cap="This is a plot of the func x2", echo=false
#theme(:bright)
#plot(f, xlim=(0,2), frame=:box, dpi=300)

fig = PyPlot.figure()
t = LinRange(-3,3,5000)
y = f.(t)

PyPlot.plot(t, y, label = "f")

PyPlot.legend()
PyPlot.xlabel("t")
PyPlot.ylabel("f")
PyPlot.title("test plot")
display(fig)

```
