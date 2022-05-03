using LinearAlgebra


import GISSMOReader

include("../src/NMRSpectraSimulator.jl")
import .NMRSpectraSimulator

import Kronecker
import Graphs
import JSON, JSON3

import Random
Random.seed!(25)

include("./helpers/operators.jl")
include("./helpers/check_mag_eq.jl")


H_params_path = "/home/roy/Documents/repo/NMRData/input/coupling_info"
dict_compound_to_filename = JSON.parsefile("/home/roy/Documents/repo/NMRData/input/compound_mapping/select_compounds.json")

#name = "L-Histidine"

unique_cs_tol = 1e-6

# load.
# name = "D-(+)-Glucose"
# load_path = joinpath(H_params_path, dict_compound_to_filename[name]["file name"])

# #valine.
# # http://gissmo.nmrfam.wisc.edu/entry/bmse000860/simulation_1
# load_path = "/home/roy/Documents/repo/NMRData/input/coupling_info/bmse000860_simulation_1.json"


# ## 1-(3,4-Difluorophenyl)ethan-1-one
# # http://gissmo.nmrfam.wisc.edu/entry/Maybridge_Ro3_Fragment_10_A08/simulation_1
# load_path = "/home/roy/Documents/repo/NMRData/input/coupling_info/Maybridge_Ro3_Fragment_10_A08_simulation_1.json"

## 2-(4-Chlorobenzylthio)-1,4,5,6-tetrahydropyrimidine hydrochloride
# http://gissmo.nmrfam.wisc.edu/entry/Maybridge_Ro3_Fragment_12_G10/simulation_1
load_path = "/home/roy/Documents/repo/NMRData/input/coupling_info/Maybridge_Ro3_Fragment_12_G10_simulation_1.json" 

# ## L-Leucine.
# # http://gissmo.nmrfam.wisc.edu/entry/bmse000042/simulation_1
# load_path = "/home/roy/Documents/repo/NMRData/input/coupling_info/bmse000042_simulation_1.json"

# ## L-isoleucine.
# # http://gissmo.nmrfam.wisc.edu/entry/bmse000041/simulation_1
# load_path = "/home/roy/Documents/repo/NMRData/input/coupling_info/bmse000041_simulation_1.json"

# ## ethanol.
# # http://gissmo.nmrfam.wisc.edu/entry/bmse000297/simulation_1
# load_path = "/home/roy/Documents/repo/NMRData/input/coupling_info/bmse000297_simulation_1.json" 


# ## D-Carnitine
# # http://gissmo.nmrfam.wisc.edu/entry/bmse000949/simulation_1
# load_path = "/home/roy/Documents/repo/NMRData/input/coupling_info/bmse000949_simulation_1.json"




H_IDs, H_css, J_IDs, J_vals = NMRSpectraSimulator.loadcouplinginfojson(load_path)


# process.
#H_inds = collect(1:length(H_IDs))
#J_inds = convertJIDstoJinds(J_IDs, H_IDs)
#g = constructspinsystemsgraphforcompound(length(H_IDs), J_inds)

J_inds_sys, J_inds_sys_local, J_IDs_sys, J_vals_sys, H_inds_sys,
    cs_sys, H_inds_singlets, cs_singlets, H_inds, J_inds,
    g = NMRSpectraSimulator.setupcsJ(H_IDs, H_css, J_IDs, J_vals)

dict_H_IDs_to_css = Dict(H_IDs .=> H_css)
dict_H_inds_to_css = Dict(H_inds .=> H_css)
dict_H_ID_to_ind = Dict(H_IDs .=> collect(1:length(H_IDs)))
dict_ind_to_H_ID = Dict( collect(1:length(H_IDs)) .=> H_IDs)

dict_J_ID_to_val = Dict(J_IDs .=> J_vals)
dict_J_ind_to_val = Dict(J_inds .=> J_vals)

dict_J_ID_to_ind = Dict(J_IDs .=> J_inds)
dict_J_ind_to_ID = Dict(J_inds .=> J_IDs)


systems_g = Graphs.connected_components(g) # node inds for spin systems.
systems_g_ID = mapQtoHID(systems_g, H_IDs)

# for manual checking.
H_IDs_sys = collect( collect( dict_ind_to_H_ID[H_inds_sys[i][k]] for k = 1:length(H_inds_sys[i])) for i = 1:length(H_inds_sys) )
H_IDs_singlets = collect( collect( dict_ind_to_H_ID[H_inds_singlets[i][k]] for k = 1:length(H_inds_singlets[i])) for i = 1:length(H_inds_singlets) )

println("H_IDs_sys: ")
display(H_IDs_sys)
println()

println("cs_sys: ")
display(cs_sys)
println()

println("H_IDs_singlets: ")
display(H_IDs_singlets)
println()

println("cs_singlets: ")
display(cs_singlets)
println()


mag_eq_sys_inds_local, mag_eq_sys_IDs,
mag_eq_sys_inds_global = getmageqcompound(g, 
H_inds_sys, dict_ind_to_H_ID, dict_H_inds_to_css,
dict_H_IDs_to_css; atol = unique_cs_tol)

println("mag_eq_sys_IDs: ")
display(mag_eq_sys_IDs)
println()

println("mag_eq_sys_inds_local: ")
display(mag_eq_sys_inds_local)
println()


#@assert 1==2

C_g = Graphs.maximal_cliques(g)

N_spin_systems = length(H_inds_sys)

i = 2

C = NMRSpectraSimulator.keeptargetintegers(C_g, H_inds_sys[i])

mag_eq_IDs0, mag_eq_inds0 = getmageqIDs(C,
    dict_ind_to_H_ID,
    dict_H_inds_to_css;
    atol = unique_cs_tol)

#
println("mag_eq_IDs0: ")
display(mag_eq_IDs0)
println()

println("mag_eq_inds0: ")
display(mag_eq_inds0)
println()

local_inds = getmageqlocalinds(H_inds_sys[i], mag_eq_inds0)

println("local_inds: ")
display(local_inds)
println()


H = combinetransitiveeqgroups(mag_eq_IDs0)
