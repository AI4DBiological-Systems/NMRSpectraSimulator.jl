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

save_root_folder = "/home/roy/MEGAsync/outputs/NMR/SH_explore"
isdir(save_root_folder) || mkpath(save_root_folder)

# molecule_names = Vector{String}(undef, 0)
# #root_folders = Vector{String}(undef, 0)
# #project_names = Vector{String}(undef, 0)
#
# push!(molecule_names, "D-(+)-Glucose")
# #push!(root_folders, "/home/roy/MEGAsync/outputs/NMR/experiments/BMRB-500-100mM")
# #push!(project_names, "D-(+)-Glucose")

#compound_select = 2

##############

# get mapping from molecule names to their spin system info json files.
H_params_path = "/home/roy/Documents/repo/NMRData/input/coupling_info"
dict_compound_to_filename = JSON.parsefile("/home/roy/Documents/repo/NMRData/input/compound_mapping/select_compounds.json")

molecule_names = collect( key for (key,val) in dict_compound_to_filename)

notes = collect( "$(i) $(molecule_names[i])" for i = 1:length(molecule_names))
outfile = joinpath(save_root_folder, "compound_numbering.txt")
open(outfile, "w") do f
  for i in notes
    println(f, i)
  end
end

#molecule_names = ["L-Serine"; "Ethanol"; "L-Isoleucine"; "L-Leucine"]
#molecule_names = ["Creatinine"; "Ethanol"]

# # machine settings taken from the BMRB 20mM 700MHz 1D 1H NOSEY experiment for glucose.
# fs = 14005.602240896402
# SW = 20.0041938620844
# ν_0ppm = 10656.011933076665

# machine and 0ppm estimate from the NRC Jan 2022 experiment.
fs = 9615.38461538462
SW = 16.0196918511501
ν_0ppm = 6753.577042707225
α_0ppm = 307479.48631657596
β_0ppm = 3.0967187803295806
λ_0ppm = 1.9474607841303089

# spin-Hamiltonian-related.
tol_coherence = 1e-2 # resonances are pairs of eigenvalues of the Hamiltonian that have quantum numbers that differ by -1. This is the tolerance away from -1 that is allowed.
α_relative_threshold = 0.05 # resonances with relative amplitude less than this factor compared to the maximum resonance in the spin group will be removed. Set to 0.0 to see every single resonance component.
Δc_partition_radius = 0.17 # determines how many resonances get grouped together. Larger number means less groups and thus more resonances per group.
SH_config_path = "/home/roy/Documents/repo/NMRData/input/SH_configs/select_compounds_SH_configs_reduce.json"

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

##############


hz2ppmfunc = uu->(uu - ν_0ppm)*SW/fs
ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)

println("Timing: mag equivalence")
@time MEs = NMRSpectraSimulator.getmageqinfomixture(molecule_names,
    H_params_path,
    dict_compound_to_filename;
    unique_cs_atol = 1e-6)

println("Timing: setupmixtureproxies()")
@time mixture_params = NMRSpectraSimulator.setupmixtureSH(molecule_names,
    H_params_path, dict_compound_to_filename, fs, SW,
    ν_0ppm;
    config_path = SH_config_path,
    MEs = MEs,
    tol_coherence = tol_coherence,
    α_relative_threshold = α_relative_threshold,
    Δc_partition_radius = Δc_partition_radius)

############## The rest are for plotting only.

#As = mixture_params[compound_select:compound_select]
As = mixture_params
dummy_SSFID = NMRSpectraSimulator.SpinSysParamsType1(0.0)

## frequency locations. For plotting.
ΩS_ppm = NMRSpectraSimulator.getPsnospininfo(As, hz2ppmfunc)
ΩS_ppm_sorted = sort(NMRSpectraSimulator.combinevectors(ΩS_ppm))

u_offset = 0.5 # in units ppm.
u_min = ppm2hzfunc(ΩS_ppm_sorted[1] - u_offset)
u_max = ppm2hzfunc(ΩS_ppm_sorted[end] + u_offset)

# This is the frequency range that we shall work with.
P = LinRange(hz2ppmfunc(u_min), hz2ppmfunc(u_max), 50000)
U = ppm2hzfunc.(P)
U_rad = U .* (2*π)

# #select_ind = findfirst(xx->xx=="Ethanol", molecule_names)
# select_ind = findfirst(xx->xx=="L-Isoleucine", molecule_names)
# A = As[select_ind]
# ME = MEs[select_ind]
#
# import Destruct
# tmp = collect( NMRSpectraSimulator.createorderingfromeqinds(ME[i], A.N_spins_sys[i]) for i = 1:length(A.N_spins_sys) )
# κs_β_ordering, κs_β_DOF = Destruct.destruct(tmp)
#
# A.Δc_m_compound
# A.part_inds_compound
# c_bar = A.Δc_bar
# c0 = A.Δc_m_compound
#
# # c = NMRSpectraSimulator.reduceΔc(A.Δc_m_compound, ME, A.N_spins_sys)
# # c2 = NMRSpectraSimulator.reduceΔc(A.Δc_m_compound[1], ME[1], A.N_spins_sys[1])
# # sum(c[1]-c2)

#@assert 1==2


include("plot.jl")
