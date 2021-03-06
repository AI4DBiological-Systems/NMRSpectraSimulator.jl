# make this into test script, spline proxy vs lorentzian.

include("../src/NMRSpectraSimulator.jl")
import .NMRSpectraSimulator

using LinearAlgebra
using FFTW
import PyPlot
#import BSON
#import JLD

#import Clustering
import Statistics

import Random
Random.seed!(25)

PyPlot.close("all")
fig_num = 1

PyPlot.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])

### user inputs.

tol_coherence = 1e-2
α_relative_threshold = 0.05
λ0 = 3.4
Δcs_max = 0.2
κ_λ_lb = 0.5
κ_λ_ub = 2.5

molecule_names = ["DSS";]
#molecule_names = ["L-Serine"; "L-Phenylalanine"; "DSS";]
#molecule_names = ["D-(+)-Glucose"; "DSS"]

# # # get the all GISSMO entries.
# import GISSMOReader
# tmp = GISSMOReader.getGISSMOentriesall()
# GISSMO_entries = GISSMOReader.extractfields(tmp, "entry")
# molecule_names = GISSMOReader.extractfields(tmp, "molecule_name")


# machine values taken from the BMRB 700 MHz 20 mM glucose experiment.
fs = 14005.602240896402
SW = 20.0041938620844
ν_0ppm = 10656.011933076665

# # machine values for the BMRB 500 MHz glucose experiment.
# ν_0ppm = 6752.490995937095
# SW = 16.0196917451925
# fs = 9615.38461538462

# path to the GISSMO Julia storage folder.
H_params_path = "/home/roy/Documents/repo/NMRData/input/coupling_info"

### end inputs.

hz2ppmfunc = uu->(uu - ν_0ppm)*SW/fs
ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)


Δcs_max_mixture = collect( Δcs_max for i = 1:length(molecule_names))

println("Timing: setupmixtureproxies()")
dummy_SSFID = NMRSpectraSimulator.SpinSysFIDType2(0.0)
@time mixture_params = NMRSpectraSimulator.setupmixtureproxies(molecule_names,
    H_params_path, ppm2hzfunc, fs, SW, λ0,
    ν_0ppm, dummy_SSFID;
    tol_coherence = tol_coherence,
    α_relative_threshold = α_relative_threshold)
As = mixture_params


#@assert 1==43

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

# purposely distort the spectra by setting non-autophased FID values.
compound_select = 1
spin_sys_select = 1
Ag = As[compound_select]
Ag.ss_params.d[spin_sys_select] = rand(length(Ag.part_inds_compound[spin_sys_select]))
Ag.ss_params.κs_λ[spin_sys_select] = rand(length(Ag.part_inds_compound[spin_sys_select])) .+ 1
Ag.ss_params.κs_β = collect( rand(length(Ag.ss_params.κs_β[i])) .* (2*π) for i = 1:length(Ag.ss_params.κs_β) )


f = uu->NMRSpectraSimulator.evalmixture(uu, mixture_params)

# test params.
ΩS_ppm = collect( hz2ppmfunc.( NMRSpectraSimulator.combinevectors(A.Ωs) ./ (2*π) ) for A in mixture_params )
#ΩS_ppm_flat = NMRSpectraSimulator.combinevectors(ΩS_ppm)




P = LinRange(hz2ppmfunc(u_min), hz2ppmfunc(u_max), 50000)
U = ppm2hzfunc.(P)
U_rad = U .* (2*π)

## parameters that affect qs.
# A.d, A.κs_λ, A.κs_β
# A.d_singlets, A.αs_singlets, A.Ωs_singlets, A.β_singlets, A.λ0, A.κs_λ_singlets
q = uu->NMRSpectraSimulator.evalitpproxymixture(uu, mixture_params)

f_U = f.(U_rad)
q_U = q.(U_rad)

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
