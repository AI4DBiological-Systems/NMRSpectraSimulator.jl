########## compare the partions between experiments at different machine settings of the same molecule.

#import NMRDataSetup

#import NMRSpectraSimulator

include("../src/NMRSpectraSimulator.jl")
import .NMRSpectraSimulator

using LinearAlgebra
using FFTW
import PyPlot
import BSON
import JLD

import Statistics
import Random


PyPlot.close("all")
fig_num = 1

Random.seed!(25)
PyPlot.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])


# helper functions.
function getΔcinfo!(part_dic, mixture_params_dic, project_name; verbose_flag = true)

    As, _ = mixture_params_dic["$(project_name)"]

    A = As[1]

    if verbose_flag
        println("$(project_name): Partition sizes:")
        tmp = collect( length(A.part_inds_compound[k]) for k = 1:length(A.part_inds_compound) )
        display(tmp)
        println("Δc_partition_radius = ", Δc_partition_radius)
        println()
    end

    for i = 1:length(A.Ωs)

        inds_set = A.part_inds_compound[i]
        cs = A.Δc_m_compound[i]
        ps = collect( cs[inds_set[i]] for i = 1:length(inds_set) )


        #centroid_ps = collect( Statistics.mean(ps[k]) for k = 1:length(ps) )
        centroid_ps = A.Δc_bar[i]

        ds = collect( collect(ps[k][l] - centroid_ps[k] for l = 1:length(ps[k])) for k = 1:length(ps))
        ns = collect( norm.(ds[k]) for k = 1:length(ds) )
        max_ns = collect( maximum(ns[k]) for k = 1:length(ns))

        part_dic["$(project_name)_spin_group_$(i)"] = Dict("max_ns" => max_ns, "centroid_ps" => centroid_ps)

        if verbose_flag
            println("Partition info for $(project_name), spin group $(i):")

            println("maximum l-2 norm distance between centroid Δc_bar and each Δc_{l} for $(project_name), spin group $(i):")

            max_ns_rounded = collect( round.(max_ns[k], sigdigits = 2) for k = 1:length(max_ns))
            display(max_ns_rounded)
            println()

            centroid_ps_rounded = collect( round.(centroid_ps[k], sigdigits = 2) for k = 1:length(centroid_ps))
            println("centroid Δc_bar for each resonance component group in this partition, rounded to 2 significant digits:")
            display(centroid_ps_rounded)

            println()
        end
    end

end

function getsimulation(projects_dir, project_name, molecule_names, tol_coherence::T, α_relative_threshold, Δc_partition_radius,
    dummy_SSFID;
    Δcs_max = 0.2, κ_λ_lb = 0.5, κ_λ_ub = 2.5) where T

    println("Now on $(project_name)")

    ##### load the metadata of the machine settings of an actual experiment.
    load_folder_path = joinpath(projects_dir, project_name)
    load_path = joinpath(load_folder_path, "experiment.bson")
    dic = BSON.load(load_path)
    fs = dic[:fs]
    SW = dic[:SW]
    ν_0ppm = dic[:ν_0ppm]
    λ_0ppm = dic[:λ_0ppm]

    hz2ppmfunc = uu->(uu - ν_0ppm)*SW/fs
    ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)
    ### end load metadata.



    ####### compute surrogate.
    Δcs_max_mixture = collect( Δcs_max for i = 1:length(molecule_names))

    # get a surrogate where K_{n,i} is encouraged to be no larger than `early_exit_part_size`.
    mixture_params = NMRSpectraSimulator.setupmixtureproxies(molecule_names,
        H_params_path, ppm2hzfunc, fs, SW,
        λ_0ppm, ν_0ppm, dummy_SSFID;
        tol_coherence = tol_coherence,
        α_relative_threshold = α_relative_threshold,
        Δc_partition_radius = Δc_partition_radius)

    # # frequency locations. For plotting.
    # ΩS_ppm = NMRSpectraSimulator.getPsnospininfo(mixture_params, hz2ppmfunc)
    # ΩS_ppm_sorted = sort(NMRSpectraSimulator.combinevectors(ΩS_ppm))
    #
    # u_offset = 0.5 # in units ppm.
    # u_min = ppm2hzfunc(ΩS_ppm_sorted[1] - u_offset)
    # u_max = ppm2hzfunc(ΩS_ppm_sorted[end] + u_offset)
    #
    # P = LinRange(hz2ppmfunc(u_min), hz2ppmfunc(u_max), 50000)
    # U = ppm2hzfunc.(P)

    return mixture_params, dic
end

function runbatchsimulations(projects_dir,
    molecule_names, tol_coherence, α_relative_threshold, Δc_partition_radius, dummy_SSFID,
    Δcs_max, κ_λ_lb, κ_λ_ub)

    mixture_params_dic = Dict()
    partition_dic = Dict()

    for project_name in project_names
        mixture_params_dic["$(project_name)"] = getsimulation(projects_dir, project_name,
        molecule_names, tol_coherence, α_relative_threshold, Δc_partition_radius, dummy_SSFID;
        Δcs_max = Δcs_max, κ_λ_lb = κ_λ_lb, κ_λ_ub = κ_λ_ub)

        getΔcinfo!(partition_dic, mixture_params_dic, project_name; verbose_flag = true)
    end

    return mixture_params_dic, partition_dic
end

##### global constants.
projects_dir = "/home/roy/MEGAsync/outputs/NMR/experiments"
H_params_path = "/home/roy/Documents/repo/NMRData/input/coupling_info"
cs_config_path = "/home/roy/Documents/repo/NMRData/input/reduced_cs_config.txt"
#####


##### user input for surrogate.

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

dummy_SSFID = NMRSpectraSimulator.SpinSysParamsType1(0.0) # level 2 model.
##### end inputs for surrogate.




##### user input for loading data.
molecule_names = ["D-(+)-Glucose";]

project_names = Vector{String}(undef, 0)

# 500 MHz.
push!(project_names, "BMRB-500-2mM-D-(+)-Glucose")
push!(project_names, "BMRB-500-0.5mM-D-(+)-Glucose")

# 600 MHZ.
push!(project_names, "BMRB-600-100mM-D-(+)-Glucose")
push!(project_names, "NRC-glucose-2018")
push!(project_names, "NRC-4_amino_acid-Jan2022-1")

# 700 MHZ.
push!(project_names, "BMRB-700-20mM-D-(+)-Glucose")
##### end inputs for loading data.

println("Working on ", molecule_names[1])
runbatchsimulations(projects_dir,
molecule_names, tol_coherence, α_relative_threshold, Δc_partition_radius, dummy_SSFID,
Δcs_max, κ_λ_lb, κ_λ_ub)
println()
println()

##### user input for loading data.
molecule_names = ["L-Serine";]

project_names = Vector{String}(undef, 0)

# 500 MHz.
push!(project_names, "BMRB-500-2mM-L-Serine")
push!(project_names, "BMRB-500-0.5mM-L-Serine")

# 600 MHZ.
push!(project_names, "NRC-4_amino_acid-Jan2022-1")

# 700 MHZ.
push!(project_names, "BMRB-700-20mM-L-Serine")
##### end inputs for loading data.

println("Working on ", molecule_names[1])
runbatchsimulations(projects_dir,
molecule_names, tol_coherence, α_relative_threshold, Δc_partition_radius, dummy_SSFID,
Δcs_max, κ_λ_lb, κ_λ_ub)
println()
println()


# @assert 1==4

# # full tutorial on accessing Δc, Δc_bar in another file. with plots.

# lens = collect( length(ps[i]) for i = 1:length(ps))

# centroid_ps_rounded = collect( round.(centroid_ps[k], digits = 2) for k = 1:length(centroid_ps))

# println("ps lengths: ", lens)
# println()


# ###### check yourself.
# α_tol = 0.01
# r = 0.1

# i_select = 1
# part_inds = NMRSpectraSimulator.partitionresonancesbyneighbors(A.Δc_m_compound[i_select],
#             A.αs[i_select], α_tol; radius = r)

# norm(A.part_inds_compound[i_select] - part_inds)
