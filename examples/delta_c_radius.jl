
# I am here. turn this into a tutorial.
# then, move onto NMRCalibrate's fit.jl, to fit n,i.
# repeat for n,i,k.
# repeat for different compounds and regions.


using LinearAlgebra, FFTW
import BSON
#import Statistics

#import NMRSpectraSimulator

include("../src/NMRSpectraSimulator.jl")
import .NMRSpectraSimulator

# for loading something with Interpolations.jl
import OffsetArrays
import Interpolations

import PlotlyJS
import Plots
Plots.plotly()

#import Destruct

#include("./helpers/final_helpers.jl")
#include("./helpers/resonance_helpers.jl")

#import PyPlot
#PyPlot.close("all")
#fig_num = 1
#PyPlot.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])

import Random
Random.seed!(25)

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

    ## uncomment for markers.
    # Plots.plot!(plot_obj, P, f.(q_U),
    # markershape = :circle,
    # seriestype = :scatter,
    # xflip = true)

    return plot_obj, q_U, qs_U, q_singlets_U
end


##### user inputs.

# simulation compounds. # paths.
plot_title = "Number of groups vs. Δc deviation tolerance"

# molecule_names = ["L-Histidine";]
# root_folder = "/home/roy/MEGAsync/outputs/NMR/experiments/BMRB-500-0.5mM"
# project_path = joinpath(root_folder, "L-Histidine")


molecule_names = ["D-(+)-Glucose";]
root_folder = "/home/roy/MEGAsync/outputs/NMR/experiments/BMRB-500-0.5mM"
project_path = joinpath(root_folder, "D-(+)-Glucose")


H_params_path = "/home/roy/Documents/repo/NMRData/input/coupling_info"

load_path = joinpath(project_path, "experiment.bson")
save_folder = joinpath(project_path, "plots")

# proxy-related.
tol_coherence = 1e-2
α_relative_threshold = 0.05
#α_relative_threshold = 0.01
#α_relative_threshold = 0.0

#Δc_partition_radius = 0.9 #1e-1
Δc_partition_radius = 0.17
#Δc_partition_radius = 1e-1

Δcs_max = 0.2 # for proxy.
κ_λ_lb = 0.5
κ_λ_ub = 2.5

dummy_SSFID = NMRSpectraSimulator.SpinSysFIDType1(0.0) # level 2 model.

# # plotting related options.
# display_flag = false
# prune_low_signal_for_display_flag = false
# display_reduction_factor = 100
# display_threshold_factor =  α_relative_threshold/10
# canvas_size = (1600, 900)
# reset_plot_dir_flag = true
# plot_imag_and_mag_flag = false
# save_plot_flag = true

##### end user inputs.

plot_title
molecule_names,
base_path_JLD
project_path
tol_coherence
α_relative_threshold,
Δc_partition_radius
Δcs_max
κ_λ_lb
κ_λ_ub
display_flag = false
prune_low_signal_for_display_flag = false
display_reduction_factor = 100
display_threshold_factor =  α_relative_threshold/10
canvas_size = (1600, 900)
reset_plot_dir_flag = true
plot_imag_and_mag_flag = false
save_plot_flag = true

load_path = joinpath(project_path, "experiment.bson")
save_folder = joinpath(project_path, "plots")
isdir(save_folder) || mkpath(save_folder)

### load block.
dic = BSON.load(load_path)
fs = dic[:fs]
SW = dic[:SW]
ν_0ppm = dic[:ν_0ppm]
λ_0ppm = dic[:λ_0ppm]

hz2ppmfunc = uu->(uu - ν_0ppm)*SW/fs
ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)


Δcs_max_mixture = collect( Δcs_max for i = 1:length(molecule_names))

##### continue here later. first, revamp the GISSMO repo, and how we call up JLD's from cpound names.
file_name = "$(molecule_name)_$(record_name)"

css_sys, J_IDs_sys, J_vals_sys, cs_singlets,
J_lb_sys, J_ub_sys, cs_lb_sys, cs_ub_sys,
cs_LUT, reference_concentration, solute_concentration,
_ = GISSMOReader.loadGISSMOmolecule(base_path, file_name)

@assert 1==2

function trydiffΔcradius(Δc_partition_radius_candidates::Vector{T},
    molecule_names, base_path_JLD, hz2ppmfunc, ppm2hzfunc,
    fs, SW, λ0, ν_0ppm, early_exit_part_size, Δcs_max, tol_coherence,
    α_relative_threshold,
    dummy_SSFID::SST)::Tuple{Vector{NMRSpectraSimulator.CompoundFIDType{T,SST}}, T} where {T <: Real, SST}

    @assert early_exit_part_size > 0
    @assert all(Δc_partition_radius_candidates .> zero(T))

    Δcs_max_mixture = collect( Δcs_max for i = 1:length(molecule_names))

    for Δc_partition_radius in Δc_partition_radius_candidates[1:end-1]

        mixture_params = NMRSpectraSimulator.setupmixtureproxies(molecule_names,
            base_path_JLD, ppm2hzfunc, fs, SW,
            λ0, ν_0ppm, dummy_SSFID;
            tol_coherence = tol_coherence,
            α_relative_threshold = α_relative_threshold,
            Δc_partition_radius = Δc_partition_radius)
        As = mixture_params

        if all( all(NMRCalibrate.displaypartitionsizes(As[n]) .<= early_exit_part_size) for n = 1:length(As) )
            return mixture_params, Δc_partition_radius
        end
    end

    mixture_params = NMRSpectraSimulator.setupmixtureproxies(molecule_names,
    base_path_JLD, ppm2hzfunc, fs, SW, λ0,
    ν_0ppm, dummy_SSFID;
    tol_coherence = tol_coherence,
    α_relative_threshold = α_relative_threshold,
    Δc_partition_radius = Δc_partition_radius_candidates[end])

    return mixture_params, Δc_partition_radius_candidates[end]
end

# get a surrogate where K_{n,i} is encouraged to be no larger than `early_exit_part_size`.
# println("Timing: setupmixtureproxies()")
# @time mixture_params = NMRSpectraSimulator.setupmixtureproxies(molecule_names,
mixture_params = NMRSpectraSimulator.setupmixtureproxies(molecule_names,
    base_path_JLD, ppm2hzfunc, fs, SW,
    λ_0ppm, ν_0ppm, dummy_SSFID;
    tol_coherence = tol_coherence,
    α_relative_threshold = α_relative_threshold,
    Δc_partition_radius = Δc_partition_radius)
As = mixture_params
A = As[1]



## frequency locations. For plotting.
ΩS_ppm = NMRSpectraSimulator.getPsnospininfo(mixture_params, hz2ppmfunc)
ΩS_ppm_sorted = sort(NMRSpectraSimulator.combinevectors(ΩS_ppm))

u_offset = 0.5 # in units ppm.
u_min = ppm2hzfunc(ΩS_ppm_sorted[1] - u_offset)
u_max = ppm2hzfunc(ΩS_ppm_sorted[end] + u_offset)

P = LinRange(hz2ppmfunc(u_min), hz2ppmfunc(u_max), 50000)
U = ppm2hzfunc.(P)
U_rad = U .* (2*π)

## end frequency locations.

# println("Timing: fitproxies()")
# @time NMRSpectraSimulator.fitproxies!(As;
NMRSpectraSimulator.fitproxies!(As;
    κ_λ_lb = κ_λ_lb,
    κ_λ_ub = κ_λ_ub,
    u_min = u_min,
    u_max = u_max,
    Δr = 1.0,
    Δκ_λ = 0.05)

# store fitted surrogate for future use.
save_simulation_path = joinpath(project_path, "simulation.bson")

#NMRSpectraSimulator.removeauxinfo!(mixture_params)
#BSON.bson(save_simulation_path, mixture_params = mixture_params)

# create the functions for each resonance group.
qs = collect( collect( ωω->A.qs[i][k](ωω-A.ss_params.d[i], A.ss_params.κs_λ[i]) for k = 1:length(A.qs[i]) ) for i = 1:length(A.qs) )
q_singlets = ωω->NMRSpectraSimulator.evalsinglets(ωω, A.d_singlets, A.αs_singlets, A.Ωs_singlets, A.β_singlets, A.λ0, A.κs_λ_singlets)

# create the function for the entire compound.
q = uu->NMRSpectraSimulator.evalitpproxymixture(uu, As[1:1])

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

#println("sanity check. This should be numerically zero: ", discrepancy)
@assert discrepancy < 1e-14

@assert 1==23

# plot.
P_display = P
U_display = U
if prune_low_signal_for_display_flag

    inds, _ = NMRSpectraSimulator.prunelowsignalentries(q_U, display_threshold_factor, display_reduction_factor)
    P_display = P[inds]
    U_display = U[inds]
end

#println("length(P) = ", length(P))
#println("length(P_display) = ", length(P_display))

# should be clear the plotting directory before saving plots?
if reset_plot_dir_flag && save_plot_flag
    paths_to_delete = readdir(save_folder, join = true)
    rm.(paths_to_delete)
end

plots_save_path = joinpath(save_folder, "groups_real.html")
title_string = "$(plot_title), $(molecule_names), real spectrum"
plot_obj, q_U, qs_U, q_singlets_U = plotgroups(title_string, P_display, U_display, q, qs, q_singlets, real, P[1]; canvas_size = canvas_size)
if save_plot_flag
    Plots.savefig(plot_obj, plots_save_path)
end
if display_flag
    display(plot_obj)
end

#TODO sample increments for Δc_partition_radius, plot its corresponding size for each partition.
# then option to use different Δc_partition_radius for each compound, spin group.
# then read in config file for each copmpound. no config entry, use default Δc_partition_radius
