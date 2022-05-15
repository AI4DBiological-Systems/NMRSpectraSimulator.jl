include("../src/NMRSpectraSimulator.jl")
import .NMRSpectraSimulator

using LinearAlgebra
using FFTW

import PlotlyJS
using Plots; plotly()

import OffsetArrays
import Interpolations

import BSON
import JSON

import Random
Random.seed!(25)


# simulation compounds. # paths.
plot_title = "resonance groupings"

# name of molecules to simulate.
# Limit to only one compound in this tutorial since we just want to visualize the different resonance groups for one compound.
#molecule_names = ["L-Histidine";]


molecule_names = ["Ethanol";]

# get mapping from molecule names to their spin system info json files.
H_params_path = "/home/roy/Documents/repo/NMRData/input/coupling_info"
dict_compound_to_filename = JSON.parsefile("/home/roy/Documents/repo/NMRData/input/compound_mapping/select_compounds.json")

# where the bson file is located.
root_folder = "/home/roy/MEGAsync/outputs/NMR/experiments/misc"
project_path = joinpath(root_folder, "bmse000297_ethanol")
load_path = joinpath(project_path, "experiment.bson")

# where to save the resultant plot. Warning: This script will create the following path if it doesn't it exist.
save_folder = joinpath(project_path, "plots")
isdir(save_folder) || mkpath(save_folder)

# spin-Hamiltonian-related.
tol_coherence = 1e-2 # resonances are pairs of eigenvalues of the Hamiltonian that have quantum numbers that differ by -1. This is the tolerance away from -1 that is allowed.
α_relative_threshold = 0.05 # resonances with relative amplitude less than this factor compared to the maximum resonance in the spin group will be removed. Set to 0.0 to see every single resonance component.
Δc_partition_radius = 0.17 # determines how many resonances get grouped together. Larger number means less groups and thus more resonances per group.
SH_config_path = "/home/roy/Documents/repo/NMRData/input/SH_configs/select_compounds_SH_configs.json"

# surrogate-related.
# default values. The alternative is to load from a config file.
Δr_default = 1.0 # the samples used to build the surrogate is taken every `Δr` radian on the frequency axis. Decrease for improved accuracy at the expense of computation resources.
Δκ_λ_default = 0.05 # the samples used to build thes urrogate for κ_λ are taken at this sampling spacing. Decrease for improved accuracy at the expense of computation resources.
Δcs_max_scalar_default = 0.2 # In units of ppm. interpolation border that is added to the lowest and highest resonance frequency component of the mixture being simulated.
κ_λ_lb_default = 0.5 # interpolation lower limit for κ_λ.
κ_λ_ub_default = 2.5 # interpolation upper limit for κ_λ.
λ0 = 3.4
surrogate_config_path = "/home/roy/Documents/repo/NMRData/input/surrogate_configs/select_compounds_surrogate_configs.json"

dummy_SSFID = NMRSpectraSimulator.SpinSysParamsType1(0.0) # level 2 model.

##### end user inputs.


### load block.
dic = BSON.load(load_path)
fs = dic[:fs]
SW = dic[:SW]
ν_0ppm = dic[:ν_0ppm]
λ_0ppm = dic[:λ_0ppm]

#ν_0ppm += 42e6

hz2ppmfunc = uu->(uu - ν_0ppm)*SW/fs
ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)

# css_sys = [[1.1696, 1.1696, 1.1696, 3.64079, 3.64079]]
# J_vals_sys = [[-12.5, -12.5, 7.094311, 7.094311, -12.5, 7.094311, 7.094311, 7.094311, 7.094311, -12.4]]
# J_inds_sys = [[(1, 2), (1, 3), (1, 4), (1, 5), (2, 3), (2, 4), (2, 5), (3, 4), (3, 5), (4, 5)]]
# N_spins_sys = [5;]
#
# intermediates_sys = NMRSpectraSimulator.prepcouplingalgorithm(N_spins_sys)
#
# i = 1
# αs, Ωs, coherence_mat_sys, eig_vals_sys, Q_sys, coherence_state_pairs_sys,
#        H_sys, H1, H2, M_sys = NMRSpectraSimulator.evalSCalgorithm!(css_sys[i],
#        J_vals_sys[i], J_inds_sys[i],
#        intermediates_sys[i],
#        ppm2hzfunc;
#        tol_coherence = 0.1)

name = molecule_names[1]
base_path = H_params_path

config_path = SH_config_path

load_path = joinpath(base_path, dict_compound_to_filename[name]["file name"])
H_IDs, H_css, J_IDs, J_vals = NMRSpectraSimulator.loadcouplinginfojson(load_path)

db_dict = JSON.parsefile(config_path)
dict = db_dict[name]

tol_coherence = dict["coherence tolerance"]
α_relative_threshold = dict["relative amplitude threshold"]
Δc_partition_radius = dict["maximum Δc deviation"]
simple_coherence_atol = dict["simple coherence absolute tolerance"]

J_inds_sys, J_inds_sys_local, J_IDs_sys, J_vals_sys, H_inds_sys,
    cs_sys, H_inds_singlets, cs_singlets, H_inds, J_inds,
    g = NMRSpectraSimulator.setupcsJ(H_IDs, H_css, J_IDs, J_vals)

N_spins_singlet = length.(H_inds_singlets)

N_spins_sys = collect( length(cs_sys[m]) for m = 1:length(cs_sys) )
intermediates_sys = NMRSpectraSimulator.prepcouplingalgorithm(N_spins_sys)

# css_sys becomes cs_sys
# cs_singlets_compact becomes cs_singlets.

αs_inp, Ωs_inp, coherence_mat_sys, eig_vals_sys, Q_sys,
coherence_state_pairs_sys, H_sys, states_sys, ms_sys,
M_sys = NMRSpectraSimulator.setupcompoundspectrum!(cs_sys,
    J_vals_sys, J_inds_sys_local, intermediates_sys,
    ppm2hzfunc, cs_singlets, N_spins_singlet, fs, SW;
    tol_coherence = tol_coherence)
#
T = Float64
αs_spin_sys = copy(αs_inp)
Ωs_spin_sys = copy(Ωs_inp)
αs_singlets = Vector{T}(undef, 0)
Ωs_singlets = Vector{T}(undef, 0)

αs, Ωs, part_inds_compound,
Δc_m_compound, Δc_bar = NMRSpectraSimulator.partitionresonances(coherence_state_pairs_sys,
ms_sys, αs_spin_sys, Ωs_spin_sys;
α_relative_threshold = α_relative_threshold,
Δc_partition_radius = Δc_partition_radius,
simple_coherence_atol = simple_coherence_atol)

M2 = collect( sum(ms_sys[1][i]) for i = 1:length(ms_sys[1]) )

cp = coherence_state_pairs_sys[1]
as = αs[1]
Fs = Ωs[1]
Ps = hz2ppmfunc.( Fs ./ (2*π) )

cb = Δc_bar[1]
cbr = collect( round.(cb[i], digits = 2) for i = 1:length(cb))
inds2 = [2; 3; 4; 13; 14; 11; 12]
inds1 = [1; 7; 8; 9; 10]

inds = part_inds_compound[1]

hcat(cbr[inds1]...) # 3 peaks.
hcat(cbr[inds2]...) # 4 peaks.

[cp[inds[1]] as[inds[1]]]

energies = collect( sum(as[inds[k]]) for k = 1:length(inds) )
[energies 1:length(energies)]

freqs = collect( collect(Ps[inds[k]]) for k = 1:length(inds) )
#[freqs 1:length(freqs)]


k_fin = findfirst(xx->isapprox(xx,0,atol = 1e-5), energies)
if typeof(k_fin) == Nothing
    k_fin = length(cbr)
end
[hcat(cbr[1:k_fin]...); collect(1:k_fin)']

[energies[1:k_fin] 1:length(energies[1:k_fin])]
#[freqs[1:k_fin] 1:length(freqs[1:k_fin])]

flags = trues(length(cbr))
flags[1:k_fin] .= false
cb2_not = cbr[flags]

tmp1 = collect( sum(cbr[k][1:3]) for k = 1:length(cbr) )
tmp2 = collect( sum(cbr[k][4:5]) for k = 1:length(cbr) )
[tmp1 tmp2 flags 1:length(cbr)]

# cp_groups = collect( cp[inds[k]] for k = 1:length(inds) )
#
# # check.
# c1 = cp_groups[1]
# m = ms_sys[1]
# Δc_m_compound[1][inds[1]] - collect( m[c1[i][1]] - m[c1[i][2]] for i = 1:length(c1) )


import Statistics
freqs = collect( collect(Ps[inds[k]]) for k = 1:length(inds) )
avg_freqs = collect( Statistics.mean(Ps[inds[k]]) for k = 1:length(inds) )
[avg_freqs energies 1:length(energies)]
println("[avg_freqs energies 1:length(energies)]")
display([avg_freqs energies 1:length(energies)])

#@assert 1==21

include("plot.jl")
