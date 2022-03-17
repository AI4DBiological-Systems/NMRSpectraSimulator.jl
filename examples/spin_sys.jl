## TODO untested with new scheme.

#include("../src/NMRSpectraSimulator.jl")
import NMRSpectraSimulator

using LinearAlgebra
using FFTW
import PyPlot
# import BSON
# import JLD

#import Clustering
import Statistics

PyPlot.close("all")
fig_num = 1

PyPlot.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])

### user inputs.

name = "DSS"
tol_coherence = 1e-2

# machine values taken from the BMRB 700 MHz 20 mM glucose experiment.
fs = 14005.602240896402
SW = 20.0041938620844
ν_0ppm = 10656.011933076665


# remove later.
base_path_JLD = "/home/roy/Documents/repo/NMRData/src/input/molecules"

### end inputs.

hz2ppmfunc = uu->(uu - ν_0ppm)*SW/fs
ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)

records = NMRSpectraSimulator.GISSMOReader.getGISSMOentriesall()
record_names = collect( records[i].molecule_name for i = 1:length(records) )
record_entries = collect( records[i].entry for i = 1:length(records) )

k = findfirst(xx->xx==name, record_names)
entries = [ records[k].entry ]
molecule_names = [records[k].molecule_name;]

record_name = entries[1]
molecule_name = molecule_names[1]

p_cs_sys, css_sys, cs_singlets, J_vals_sys, J_IDs_sys, intermediates_sys,
cs_LUT, cs_singlets_compact, N_spins_singlet, css, p_compound,
cs_compound = NMRSpectraSimulator.fetchSHparameters(molecule_name, record_name, base_path_JLD)


N_spins_compound = collect( length(css_sys[i]) for i = 1:length(css_sys))
# TODO take care of singlets, this phase scheme.
# TODO store metabolite q's.
αs_inp, Ωs_inp, coherence_mat_sys, eig_vals_sys, Q_sys,
coherence_state_pairs_sys, H_sys, states_sys, ms_sys, M_sys = NMRSpectraSimulator.setupcompoundspectrum!(css_sys,
p_cs_sys, J_vals_sys, J_IDs_sys, intermediates_sys, cs_LUT, ppm2hzfunc,
cs_singlets_compact, N_spins_singlet, fs, SW; tol_coherence = tol_coherence)



Ωs_ppm = hz2ppmfunc.( NMRSpectraSimulator.combinevectors(Ωs_inp) ./ (2*π) )
min_ppm = minimum(Ωs_ppm) - 0.5
max_ppm = maximum(Ωs_ppm) + 0.5

P = LinRange(min_ppm, max_ppm, 50000)
U = ppm2hzfunc.(P)

k = findfirst(xx->length(xx)==1, αs_inp)
αs_spin_sys = copy(αs_inp)
Ωs_spin_sys = copy(Ωs_inp)
if typeof(k) == Int
    αs_spin_sys = αs_spin_sys[1:k-1]
    Ωs_spin_sys = Ωs_spin_sys[1:k-1]
end


α_relative_threshold = 0.05
αs, Ωs, part_inds_compound,
Δc_m_compound = NMRSpectraSimulator.partitionresonances(coherence_state_pairs_sys,
    ms_sys, αs_spin_sys, Ωs_spin_sys; α_relative_threshold = α_relative_threshold)

# test params.
λ0 = 3.4
κ_λ_lb = 0.5
κ_λ_ub = 2.5

# κs_λ = ones(length(αs))
# κs_β = collect( zeros(N_spins_compound[i]) for i = 1:length(N_spins_compound))
# d = zeros(T, length(αs))

κs_λ = collect( NMRSpectraSimulator.convertcompactdomain(rand(), 0.0, 1.0, κ_λ_lb, κ_λ_ub) for i = 1:length(αs))
κs_β = collect( randn(N_spins_compound[i]) for i = 1:length(N_spins_compound))
d = rand(length(αs))


# oracle.
f = uu->NMRSpectraSimulator.evalcLcompoundviapartitions(uu, d,
αs, Ωs, κs_λ, κs_β, λ0, Δc_m_compound, part_inds_compound)

# proxy.
Δcs_max = 0.2
d_max = ppm2hzfunc(Δcs_max)-ppm2hzfunc(0.0)
@time qs, gs = NMRSpectraSimulator.setupcompoundpartitionitp(d_max,
    Δc_m_compound,
    part_inds_compound,
    αs, Ωs,
    λ0, minimum(U), maximum(U);
    κ_λ_lb = 0.5,
    κ_λ_ub = 2.5)

q = uu->NMRSpectraSimulator.evalitpproxysys(qs, uu, d, κs_λ, κs_β)



f_U = f.(U)
q_U = q.(U)



discrepancy = abs.(f_U-q_U)
max_val, ind = findmax(discrepancy)
println("relative discrepancy = ", norm(discrepancy)/norm(f_U))
println("max discrepancy: ", max_val)
println()

# TODO next, code up updatecs, β and λ.

import BenchmarkTools
BenchmarkTools.@btime f(U[ind]) # 42 microsec.
BenchmarkTools.@btime q(U[ind]) # 6 microsec.

import PyPlot
## visualize.
PyPlot.plot(P, real.(f_U), label = "f")
PyPlot.plot(P, real.(q_U), label = "q")

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("f vs q")
