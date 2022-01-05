include("../src/NMRSpectraSimulator.jl")
import .NMRSpectraSimulator

using LinearAlgebra
using FFTW
import PyPlot
import BSON
import JLD

#import Clustering
import Statistics

PyPlot.close("all")
fig_num = 1

PyPlot.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])

### user inputs.

#molecule_names = ["L-Histidine"; "L-Alanine"; "DSS"]
#molecule_names = ["L-Histidine"; "L-Alanine"]
molecule_names = ["L-Serine";
"L-Phenylalanine";
"DSS";]

tol_coherence = 1e-2
α_relative_threshold = 0.05
λ0 = 3.4
Δcs_max = 0.2
κ_λ_lb = 0.5
κ_λ_ub = 2.5

# machine values taken from the BMRB 700 MHz 20 mM glucose experiment.
# fs = 14005.602240896402
# SW = 20.0041938620844
# ν_0ppm = 10656.011933076665

ν_0ppm = 6752.490995937095
SW = 16.0196917451925
fs = 9615.38461538462


# remove later.
base_path_JLD = "/home/roy/Documents/data/NMR/NMRData/src/input/molecules"

### end inputs.

hz2ppmfunc = uu->(uu - ν_0ppm)*SW/fs
ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)


Δcs_max_mixture = collect( Δcs_max for i = 1:length(molecule_names))

mixture_params = NMRSpectraSimulator.setupmixtureproxies(molecule_names,
    base_path_JLD, Δcs_max_mixture, hz2ppmfunc, ppm2hzfunc, fs, SW, λ0,
    ν_0ppm;
    tol_coherence = tol_coherence,
    α_relative_threshold = α_relative_threshold)
As = mixture_params



u_min = ppm2hzfunc(-0.5)
u_max = ppm2hzfunc(4.0)

NMRSpectraSimulator.fitproxies!(As;
    κ_λ_lb = κ_λ_lb,
    κ_λ_ub = κ_λ_ub,
    u_min = u_min,
    u_max = u_max,
    Δr = 1.0,
    Δκ_λ = 0.05)


### plot.

f = uu->NMRSpectraSimulator.evalmixture(uu, mixture_params)

# test params.
ΩS_ppm = collect( hz2ppmfunc.( NMRSpectraSimulator.combinevectors(A.Ωs) ./ (2*π) ) for A in mixture_params )
#ΩS_ppm_flat = NMRSpectraSimulator.combinevectors(ΩS_ppm)




P = LinRange(hz2ppmfunc(u_min), hz2ppmfunc(u_max), 50000)
U = ppm2hzfunc.(P)

q = uu->NMRSpectraSimulator.evalitpproxymixture(uu, mixture_params)

f_U = f.(U)
q_U = q.(U)

discrepancy = abs.(f_U-q_U)
max_val, ind = findmax(discrepancy)
println("relative discrepancy = ", norm(discrepancy)/norm(f_U))
println("max discrepancy: ", max_val)
println()

## visualize.
PyPlot.plot(P, real.(f_U), label = "f")
PyPlot.plot(P, real.(q_U), label = "q")

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("f vs q")