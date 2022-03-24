"""
title: NMRSpectraSimulator.jl Documentation
author Roy Wang
date: 2022-03-23
"""

# Library description

Given a list of user-specified compounds and parameters, this library provides the simulated resonance components of the 1D 1H NMR experiment signal.

# Background: Signal Model
The signal model is taken to be the complex-valued free-induction decay (FID) model 

$f_{\text{FID}}(t,w,\alpha,\Omega,\lambda,\beta)=\sum_{n=1}^{N}w_{n}\sum_{l=1}^{L_{n}}\alpha_{n,l}e^{\iota\left(\Omega_{n,l}t-\beta_{n,l}\right)}e^{-\lambda_{n,l}t},\,t>0,$
and its Fourier transform, also referred to as the Fourier FID or FID spectra model 

$f_{\text{FFID}}(u,w,\alpha,\Omega,\lambda,\beta)=\sum_{n=1}^{N}w_{n}\sum_{l=1}^{L_{n}}\frac{\alpha_{n,l}e^{\iota\beta_{n,l}}}{\lambda_{n,l}+\iota(2\pi u-\Omega_{n,l})}.$

Here, $n$ is the compound index, $\left(n,l\right)$ is the resonance
component index for compound $n$, $N$ is the number of compounds
in the isotropic mixture sample, $L_{n}$ is the number of resonance
components for compound $n$, and $\iota$ is the imaginary number.
The parameters of this signal model are: The compound molar concentration
$w$, the resonance frequency $\Omega$, the resonance amplitude $\alpha$,
the resonance T2 decay $\lambda$, the resonance phase $\beta$. The
resonance frequencies and amplitudes of a compound can be simulated
if the physical chemistry parameters, chemical shift and J-coupling,
are known. This library uses the physical chemistry parameters from the [GISSMO website](http://gissmo.nmrfam.wisc.edu/library), which were obtained by empirically fitting a different forward model than the one used in this library against real NMR experiments. For GISSMO-related details, see [Dashti et. al 2017](https://pubs.acs.org/doi/abs/10.1021/acs.analchem.8b02660) and [Dashti et. al 2018](http://pubs.acs.org/doi/abs/10.1021/acs.analchem.7b02884) for their official journal publications.

# Background: Calibration Challenge
One approach to fitting this model to experimental data, i.e., to
calibrate this model to data, is to assign bound constraints on every
model parameter. The following forward model used for fitting against
a cost function illustrates this approach: 

$g_{\text{FFID}}(u,w,d,\lambda,\beta)=\sum_{n=1}^{N}w_{n}\sum_{l=1}^{L_{n}}\frac{\alpha_{n,l}e^{\iota\beta_{n,l}}}{\lambda_{n,l}+\iota(2\pi u+d_{n,l}-\Omega_{n,l})}.$

here, the frequency variable $u$ is in Hz. Many publicly-available
NMR spectra model fitting routines are based on an internal or external
library that contain a set of resonance frequencies $\Omega$ and
amplitudes $\alpha$ for each library compound. These values are usually
empirically fitted against experiments called compound standards,
where only the compound and the bare minimum ingredients are placed
in the mixture. In a non-trivial NMR experiment, the different compounds
interact with each other, causing small frequency shifts away from
the set of library frequencies $\Omega$. Thus, the parameter $d$
is required to compensate for this shift. It is often assumed that
the set of library resonance amplitudes $\alpha$ does not change
across different experiments. The fit against data in the calibration
problem is over the variables $u$, $w$, $d$, $\lambda$, and $\beta$.
Due to the large number of resonance components for each compound,
there are too many degrees of freedom to produce a good and sensible
fit for most NMR applications in a reasonable computation time. From
physical chemistry theory, some resonance components should shift
only in a certain manner. See figures x and y in our survey for an
example of poor and non-sensible fit of select fitting algorithms
against a real-world mixture.

# NMRSPectraSimulator: Proposed Model & Surrogate Model
Let $\left[L\right]:=1,\dots,L$, which is a set of size $L$. Our
modified FID spectra model 

$g_{mFFID}(u,d,\kappa,w,\lambda,\beta)=\sum_{n=1}^{N}w_{n}\sum_{i=1}^{G_{n}}\sum_{k=1}^{K_{n,i}}\sum_{l=1}^{L_{n,i,k}}\frac{\kappa_{n,i,k}\alpha_{n,i,k,l}e^{\iota\left\langle \beta_{n,i,k},c_{n,i,k}\right\rangle }}{\lambda_{n,i,k}+\iota(2\pi u+d_{n,i}-\Omega_{n,i,k,l})}$

have grouping information that reduce the number of variables to
calibrate. It essentially replace the index $l\in\left[L_{n}\right]$
by multi-indices $\left(i,k,l\right)\in\left[G_{n}\times K_{n,i}\times L_{n,i,k}\right]$.
Note that the size of $\left[L_{n}\right]$ is the same as $\left[G_{n}\times K_{n,i}\times L_{n,i,k}\right]$
because they describe two different ways of grouping the resonance
components of compound $n$. Instead of having all parameters indexed
by $\left(i,k,l\right)$, the frequency shift parameters $d$ is now
indexed by $\left(n,i\right)$, and both the T2 decay $\lambda$ and
phase $\beta$ are now indexed by $\left(n,i,k\right)$. The physical
interpretation and derivation of this grouping and associated
variables will be discussed in detail in our upcoming journal article.


To further speed up the computation, we use a spline-based surrogate
$q$ for these models instead of the modified model $g_{\text{mFFID}}$
during numerical optimization or statistical inference. Our local
$k_{i_{n}}$ surrogate model

$q_{n,i,k}(r,h)\approx\sum_{l=1}^{L_{n,i,k}}\frac{\alpha_{n,i,k,l}}{h+\iota(r-\Omega_{n,i,k,l})}$

is comprised of spline interpolation fits of the RHS of the above
equation; one for the real part, and one for the imaginary part. The
Interpolations.jl spline library is used. The surrogate model that
we use for the entire mixture is 

$q(u,w,\kappa,\lambda,\beta)=\sum_{n=1}^{N}w_{n}\sum_{i=1}^{G_{n}}\sum_{k=1}^{K_{n,i}}z_{n,i,k}q_{n,i,k}\left(2\pi u-d_{n,i},\lambda_{n,i,k}\right)$
where 
$z_{n,i,k}:=\kappa_{n,i,k}e^{\iota\left\langle \beta_{n,i,k},c_{n,i,k}\right\rangle }\in\mathbb{C}$ 
was used to reduce clutter.

# under construction.

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
