
import PyPlot

# include("../src/NMRSpectraSimulator.jl")
# using .NMRSpectraSimulator

using Distributed, LinearAlgebra, FFTW, SharedArrays, OffsetArrays
import DelimitedFiles, Unicode, FileIO, Random, Printf, PyCall, FiniteDiff, Optim, NLopt, Interpolations, BSON, JLD, Kronecker, BenchmarkTools, NearestNeighbors, HTTP, JSON3, Dates
include("../src/SH/SH_front_end.jl")
include("../src/SH/molecule.jl")
include("../src/SH/SH.jl")
include("../src/SH/Hamiltonian.jl")
include("../src/SH/operators.jl")
include("../src/SH/fuse_components.jl")
include("../src/SH/full_solve.jl")
include("../src/IO/load_configs.jl")
include("../src/database/GISSMO_entries.jl")
include("../src/database/database_helpers.jl")
include("../src/model/parse.jl")
include("../src/misc/utilities.jl")

using FFTW
import PyPlot
import NMRData

PyPlot.close("all")
fig_num = 1

PyPlot.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])

Random.seed!(25)

#cs_config_path = "../configs/cs_config.txt"
save_folder = "../output/experiments"

base_path_JLD = NMRData.getpath("input", "molecules")
base_path_GISSMO = NMRData.getpath("input", "GISSMO_data")

target_names_path = "/home/roy/Documents/repo/NMRSpectraSimulator/configs/DMEM_targeted_metabolites.txt"
cs_config_path = "/home/roy/Documents/repo/NMRSpectraSimulator/configs/cs_config.txt"
base_path_JLD = "/home/roy/.julia/packages/NMRData/6eos7/src/input/molecules"
base_path_GISSMO = "/home/roy/.julia/packages/NMRData/6eos7/src/input/GISSMO_data"
save_name = "test_case"
save_folder ="/home/roy/MEGAsync/outputs/NMR/dev"

experiment_full_path = "/home/roy/.julia/packages/NMRData/6eos7/src/experiments/NRC/misc/dmem_medium/Oct-22-2012"



### load_experiment.jl
target_names = readtargetnames(target_names_path)

N_targets = length(target_names)
tmp = collect( getGISSMOentry(target_names[i]) for i = 1:N_targets )
entries = collect( tmp[i].entry for i = 1:N_targets )
molecule_names = collect( tmp[i].molecule_name for i = 1:N_targets )



# # load the experiment.
# import NMRMetaboliteQuantifier
# s_t, S, hz2ppmfunc, ppm2hzfunc, ν_0ppm, fs, SW,
# α_0ppm, β_0ppm, λ_0ppm, Ω_0ppm,
# results_0ppm, dic, α_solvent, β_solvent, λ_solvent, Ω_solvent,
# results_solvent = NMRMetaboliteQuantifier.loadspectrum(experiment_full_path;
# solvent_ppm = 4.8,
# solvent_window_ppm = 0.1)
fs = 7211.53846153846
SW = 12.0143685542768
ν_0ppm = 4273.112175430351

# The following formula converts between ppm (pp) to Hz (uu).
hz2ppmfunc = uu->(uu - ν_0ppm)*SW/fs
ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)

# set up SC and proxy parameters.
tol_coherence = 1e-2 # 1e-7 #for -1 coherence.
prune_a_tol_factor = 0.01
fuse_radius_δppm = 0.001
partition_radius_δppm = 0.1
prune_zero_tol = 1e-8
Minkowski_parameter = 3.5
partition_Minkowski_parameter = 3.5
add_reference_and_solvent_flag = true
solvent_ppm_guess = 4.8
solvent_window_ppm = 0.1
update_α_flag = false

fuse_radius = ppm2hzfunc(fuse_radius_δppm) - ppm2hzfunc(0.0)
partition_radius = ppm2hzfunc(partition_radius_δppm)-ppm2hzfunc(0.0)


αs_set, Ωs_set, updateαΩFSfunc_set, p0_set,
N_p_cs_sys_set, p0_group_set = setupfit1molecules(base_path_JLD,
entries,
molecule_names,
ppm2hzfunc,
fs,
SW;
tol_coherence = tol_coherence,
prune_a_tol_factor = prune_a_tol_factor,
prune_zero_tol = prune_zero_tol,
fuse_radius = fuse_radius,
Minkowski_parameter = Minkowski_parameter,
partition_radius = partition_radius,
partition_Minkowski_parameter = partition_Minkowski_parameter,
update_α_flag = update_α_flag)

@assert 1==2
NMRSpectraSimulator.MVPfrontend(experiment_full_path,
target_names_path,
cs_config_path,
base_path_JLD,
base_path_GISSMO,
save_name,
save_folder)
