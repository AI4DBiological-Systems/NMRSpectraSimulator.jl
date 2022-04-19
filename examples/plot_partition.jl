
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
    P::LinRange{T},
    U,
    q,
    qs,
    q_singlets,
    f::Function;
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

    return plot_obj, q_U, qs_U, q_singlets_U
end


##### user inputs.

# simulation compounds. # paths.
plot_title = "resonance groupings"

# molecule_names = ["L-Histidine";]
# root_folder = "/home/roy/MEGAsync/outputs/NMR/experiments/BMRB-500-0.5mM"
# project_path = joinpath(root_folder, "L-Histidine")


molecule_names = ["D-(+)-Glucose";]
root_folder = "/home/roy/MEGAsync/outputs/NMR/experiments/BMRB-500-0.5mM"
project_path = joinpath(root_folder, "D-(+)-Glucose")


# molecule_names = ["Glycine";]
# root_folder = "/home/roy/MEGAsync/outputs/NMR/experiments/BMRB-500-0.5mM"
# project_path = joinpath(root_folder, "Glycine")




base_path_JLD = "/home/roy/Documents/repo/NMRData//input/molecules"

load_path = joinpath(project_path, "experiment.bson")
save_folder = joinpath(project_path, "plots")

# proxy-related.
tol_coherence = 1e-2
α_relative_threshold = 0.05
#α_relative_threshold = 0.01
#α_relative_threshold = 0.0

#Δc_partition_radius = 0.3 #1e-1
Δc_partition_radius = 0.17
#Δc_partition_radius = 1e-1

Δcs_max = 0.2 # for proxy.
κ_λ_lb = 0.5
κ_λ_ub = 2.5

dummy_SSFID = NMRSpectraSimulator.SpinSysFIDType1(0.0) # level 2 model.

##### end user inputs.

"""
only plots the first compound in molecule_names.
"""
function plotgroupsfullscript(plot_title, molecule_names,
   base_path_JLD, load_path, save_folder, tol_coherence, α_relative_threshold,
   Δc_partition_radius, Δcs_max, κ_λ_lb, κ_λ_ub;
   display_flag = false)
        
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

    # get a surrogate where K_{n,i} is encouraged to be no larger than `early_exit_part_size`.
    #println("Timing: setupmixtureproxies()")
    #@time mixture_params = NMRSpectraSimulator.setupmixtureproxies(molecule_names,
    mixture_params = NMRSpectraSimulator.setupmixtureproxies(molecule_names,
        base_path_JLD, Δcs_max_mixture, hz2ppmfunc, ppm2hzfunc, fs, SW,
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


    # plot.
    canvas_size = (1600, 900)

    plots_save_path = joinpath(save_folder, "groups_real.html")
    title_string = "$(plot_title), real"
    plot_obj, q_U, qs_U, q_singlets_U = plotgroups(title_string, P, U, q, qs, q_singlets, real; canvas_size = canvas_size)
    Plots.savefig(plot_obj, plots_save_path)
    if display_flag
        display(plot_obj)
    end

    plots_save_path = joinpath(save_folder, "groups_imag.html")
    title_string = "$(plot_title), imag"
    plot_obj, q_U, qs_U, q_singlets_U = plotgroups(title_string, P, U, q, qs, q_singlets, imag; canvas_size = canvas_size)
    Plots.savefig(plot_obj, plots_save_path)
    if display_flag
        display(plot_obj)
    end

    plots_save_path = joinpath(save_folder, "groups_modulus.html")
    title_string = "$(plot_title), modulus"
    plot_obj, q_U, qs_U, q_singlets_U = plotgroups(title_string, P, U, q, qs, q_singlets, abs; canvas_size = canvas_size)
    Plots.savefig(plot_obj, plots_save_path)
    if display_flag
        display(plot_obj)
    end
end

#assert 1==2

"""
root_folder should contain only folders of 1D 1H NMR compound standard experiments, with the compound name being the subfolder name.
I.e., a folder containing an experiment of L-Serine and L-Histidine should have two folders in it with those names that contain the corresponding experimental data.
"""
function batchplotgroups(plot_title, root_path,
    base_path_JLD, tol_coherence,
    α_relative_threshold, Δc_partition_radius, Δcs_max, κ_λ_lb, κ_λ_ub;
    display_flag = false)

    # get a list of experiment paths and compounds names. Assume the folder name of the experiment contains the compound.
    tmp = readdir(root_path, join = true)
    inds = findall(xx->isdir(xx), tmp)

    project_paths = tmp[inds]
    experiment_names = readdir(root_path)[inds]
    project_names = collect( "$(experiment_names[i])" for i = 1:length(experiment_names))

    for i = 1:length(project_names)

        project_path = project_paths[i]
        compound_name = experiment_names[i]

        load_path = joinpath(project_path, "experiment.bson")
        save_folder = joinpath(project_path, "plots")
        println("Now on ", project_path)

        # make sure this compound is in our library.
        records = GISSMOReader.getGISSMOentriesall()
        record_names = collect( records[i].molecule_name for i = 1:length(records) )

        k = findfirst(xx->xx==compound_name, record_names)
        
        # decide if we should simulate and plot.
        if typeof(k) != Nothing
            plotgroupsfullscript(plot_title, [compound_name;],
                base_path_JLD, load_path, save_folder, tol_coherence, α_relative_threshold,
                Δc_partition_radius, Δcs_max, κ_λ_lb, κ_λ_ub;
                display_flag = display_flag)
        else
            println("Compound not in the local GISSMO library. Skip.")
        end
    end

end


import GISSMOReader

root_folder = "/home/roy/MEGAsync/outputs/NMR/experiments/BMRB-500-0.5mM"
root_folder = "/home/roy/MEGAsync/outputs/NMR/experiments/BMRB-500-2mM"
root_folder = "/home/roy/MEGAsync/outputs/NMR/experiments/BMRB-700-20mM"
batchplotgroups(plot_title, root_folder,
    base_path_JLD, tol_coherence,
    α_relative_threshold, Δc_partition_radius, Δcs_max, κ_λ_lb, κ_λ_ub)

####