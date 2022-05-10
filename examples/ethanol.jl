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
isdir(save_folder) || mkdir(save_folder)

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

ν_0ppm += 1e6

hz2ppmfunc = uu->(uu - ν_0ppm)*SW/fs
ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)


# get a surrogate where K_{n,i} is encouraged to be no larger than `early_exit_part_size`.
#println("Timing: setupmixtureproxies()")
#@time mixture_params = NMRSpectraSimulator.setupmixtureproxies(molecule_names,
mixture_params = NMRSpectraSimulator.setupmixtureSH(molecule_names,
    H_params_path, dict_compound_to_filename, fs, SW,
    ν_0ppm;
    config_path = SH_config_path,
    tol_coherence = tol_coherence,
    α_relative_threshold = α_relative_threshold,
    Δc_partition_radius = Δc_partition_radius)
As = mixture_params

# We only work with a single compound in this tutorial. Assign a new object for this compound to reduce clutter.
A = As[1];
c2 = A.Δc_m_compound[1]
c2r = collect( round.(c2[i], digits = 2) for i = 1:length(c2))
c2a = collect( abs.(c2[i]) for i = 1:length(c2))
sum_c2a = collect(sum(c2a[i]) for i = 1:length(c2a))

flags = isapprox.(sum_c2a, 1.0, atol = 1e-3)
count(flags)

c2a[flags]

b2 = A.Δc_bar[1]
b2r = collect( round.(b2[i], digits = 2) for i = 1:length(b2))

#
collect( sum(c) for c in A.Δc_m_compound[1] )

N_p_actives = collect( count(isapprox.(c, 1.0, atol = 1e-1)) for c in A.Δc_m_compound[1] )
N_n_actives = collect( count(isapprox.(c, -1.0, atol = 1e-1)) for c in A.Δc_m_compound[1] )

count(N_p_actives .> 2)

count(N_p_actives .> 0)
count(N_n_actives .> 0)

count(N_p_actives .== 1)
count(N_n_actives .== 1)


N_passive = collect( count(isapprox.(c, 0.0, atol = 2e-1)) for c in A.Δc_m_compound[1] )

count(N_passive .== 4)

[N_p_actives N_n_actives ]

## frequency locations. For plotting.
ΩS_ppm = NMRSpectraSimulator.getPsnospininfo(mixture_params, hz2ppmfunc)
ΩS_ppm_sorted = sort(NMRSpectraSimulator.combinevectors(ΩS_ppm))

u_offset = 0.5 # in units ppm.
u_min = ppm2hzfunc(ΩS_ppm_sorted[1] - u_offset)
u_max = ppm2hzfunc(ΩS_ppm_sorted[end] + u_offset)

# This is the frequency range that we shall work with.
P = LinRange(hz2ppmfunc(u_min), hz2ppmfunc(u_max), 50000)
U = ppm2hzfunc.(P)
U_rad = U .* (2*π)

# fit the surrogate.
Bs = NMRSpectraSimulator.fitproxies(As, dummy_SSFID, λ0;
    names = molecule_names,
    config_path = surrogate_config_path,
    Δcs_max_scalar_default = Δcs_max_scalar_default,
    κ_λ_lb_default = κ_λ_lb_default,
    κ_λ_ub_default = κ_λ_ub_default,
    u_min = u_min,
    u_max = u_max,
    Δr_default = Δr_default,
    Δκ_λ_default = Δκ_λ_default)
B = Bs[1] # looking at the first (and only) molecule in this example script.

# create the functions for each resonance group.
qs = collect( collect( ωω->B.qs[i][k](ωω-B.ss_params.d[i], B.ss_params.κs_λ[i]) for k = 1:length(B.qs[i]) ) for i = 1:length(B.qs) )
q_singlets = ωω->NMRSpectraSimulator.evalsinglets(ωω, B.d_singlets, A.αs_singlets, A.Ωs_singlets, B.β_singlets, B.λ0, B.κs_λ_singlets)

# create the function for the entire compound.
q = uu->NMRSpectraSimulator.evalitpproxymixture(uu, As[1:1], Bs[1:1])

# evaluate at the plotting positions.
q_U = q.(U_rad)

qs_U = collect( collect( qs[i][k].(U_rad) for k = 1:length(qs[i]) ) for i = 1:length(qs) )
q_singlets_U = q_singlets.(U_rad)


# sanity check.
q_check_U = q_singlets_U
if !isempty(qs) # some compounds only have singlets.
    q_check_U += sum( sum( qs[i][k].(U_rad) for k = 1:length(qs[i]) ) for i = 1:length(qs) )
end
discrepancy = norm(q_check_U- q_U)
println("sanity check. This should be numerically zero: ", discrepancy)


"""
save resonance groupings of a compound.
Choices for `f`` are: `real()`, `imag()`, or `abs()`. These corresponds to real part, imaginary part, and magnitude spectrum, respectively.
"""
function plotgroups(title_string::String,
    P,
    U,
    q,
    qs,
    q_singlets,
    f::Function,
    return_val_type::T;
    canvas_size::Tuple{Int,Int} = (1000, 400)) where T

    U_rad = U .* (2*π)

    q_U = q.(U_rad)
    plot_obj = Plots.plot( P,
        f.(q_U),
        title = title_string,
        label = "sum of sub-models",
        seriestype = :line,
        ticks = :native,
        xlims = (P[1],P[end]),
        hover = P,
        linewidth = 4,
        xflip = true,
        size = canvas_size)

    qs_U = Vector{Vector{Vector{Complex{T}}}}(undef, length(qs))
    for i = 1:length(qs)

        qs_U[i] = Vector{Vector{Complex{T}}}(undef, length(qs[i]))
        for k = 1:length(qs[i])

            qs_U[i][k] = qs[i][k].(U_rad)

            Plots.plot!(plot_obj, P, f.(qs_U[i][k]), label = "sub-group ($(i),$(k))",
                seriestype = :line,
                linestyle = :dot,
                xflip = true,
                linewidth = 4)
        end
    end

    q_singlets_U = q_singlets.(U_rad)
    Plots.plot!(plot_obj, P, f.(q_singlets_U), label = "group of all singlets",
        seriestype = :line,
        linestyle = :dot,
        xflip = true,
        linewidth = 4)
    #
    #Plots.plot!(plot_obj, P, f.(q_U),
    #markershape = :circle,
    #seriestype = :scatter,
    #xflip = true)

    return plot_obj, q_U, qs_U, q_singlets_U
end



# reduce the plotting positions for low signal regions. Otherwise the plot store size will be too large, and the time to load the plot will be long.
display_reduction_factor = 100
display_threshold_factor =  α_relative_threshold/10

inds, _ = NMRSpectraSimulator.prunelowsignalentries(q_U, display_threshold_factor, display_reduction_factor)
P_display = P[inds]
U_display = U[inds]

# plot.
canvas_size = (1000, 400)

plots_save_path = joinpath(save_folder, "groups_real.html")
title_string = "$(plot_title), real"
plot_obj, q_U, qs_U, q_singlets_U = plotgroups(title_string, P_display, U_display, q, qs, q_singlets, real, P[1]; canvas_size = canvas_size)
Plots.savefig(plot_obj, plots_save_path)
display(plot_obj)
