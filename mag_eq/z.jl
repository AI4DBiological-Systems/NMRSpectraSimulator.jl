using LinearAlgebra


import GISSMOReader

include("../src/NMRSpectraSimulator.jl")
import .NMRSpectraSimulator

import Kronecker
import Graphs
import JSON, JSON3

include("./helpers/operators.jl")
include("./helpers/check_mag_eq.jl")


H_params_path = "/home/roy/Documents/repo/NMRData/input/coupling_info"
dict_compound_to_filename = JSON.parsefile("/home/roy/Documents/repo/NMRData/input/compound_mapping/select_compounds.json")

#name = "L-Histidine"

unique_cs_tol = 1e-6
unique_cs_digits = 5
zero_tol_sigdigits = 6

# load.
#name = "D-(+)-Glucose"
#load_path = joinpath(H_params_path, dict_compound_to_filename[name]["file name"])

#load_path = "/home/roy/Documents/repo/NMRData/input/coupling_info/bmse000860_simulation_1.json" # http://gissmo.nmrfam.wisc.edu/entry/bmse000860/simulation_1
load_path = "/home/roy/Documents/repo/NMRData/input/coupling_info/bmse000297_simulation_1.json" # http://gissmo.nmrfam.wisc.edu/entry/bmse000297/simulation_1
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

C = Graphs.maximal_cliques(g)
C_IDs = collect( collect( dict_ind_to_H_ID[C[i][j]] for j = 1:length(C[i]) ) for i = 1:length(C) )


# ###
# T = Float64
# i = 1

# ### partition the nucleus by unique cs.
# cs = collect( dict_H_inds_to_css[C[i][k]] for k = 1:length(C[i]) )
# cs_IDs = collect( dict_ind_to_H_ID[C[i][k]] for k = 1:length(C[i]) )

# unique_cs = unique(round.(cs, digits = unique_cs_digits))
# unique_cs_inds = collect( inds = findall(xx->isapprox(unique_cs[k], xx; atol = unique_cs_tol), cs) for k = 1:length(unique_cs) )

# ### for each partition subset that contains 2 nuclei, check condition for mag eq.
# # I am here. function for determining whether given 2+ suspected eq nuclei and their interaction nuclei list.
# k = 2
# # common_cs_IDs = cs_IDs[unique_cs_inds[k]]
# # p_IDs = getpairs(common_cs_IDs)

# # ## test common J within pairs.
# # J_p = collect( getJfromdict(p_IDs[l][1], p_IDs[l][2], dict_J_ID_to_val) for l = 1:length(p_IDs) )
# # pass_common_J_flag = isallsame(J_p; atol = unique_cs_tol)

# # ## test common J with other IDs.
 
# # # find all connections to first pair.
# # t_IDs = getJIDstest(J_IDs, common_cs_IDs, dict_H_IDs_to_css; atol = unique_cs_tol)

# # # check J-values.
# # J_t = collect( getJfromdict(t_IDs[l][1], t_IDs[l][2], dict_J_ID_to_val) for l = 1:length(t_IDs) )
# # pass_test_J_flag = isallsame(J_t; atol = unique_cs_tol)

# pass_flags = collect( checkmageq(unique_cs_inds[k], cs_IDs, dict_J_ID_to_val; atol = unique_cs_tol) for k = 1:length(unique_cs_inds) )
# common_cs_IDs = collect( cs_IDs[unique_cs_inds[k]] for k = 1:length(unique_cs_inds) )
# Q = common_cs_IDs[pass_flags]

C_IDs_mag_eq = getmageqIDs(C,
    dict_ind_to_H_ID,
    dict_H_inds_to_css;
    cs_round_digits = unique_cs_digits,
    atol = unique_cs_tol)

println("C_IDs_mag_eq: ")
display(C_IDs_mag_eq)
println()

println("C_IDs: ")
display(C_IDs)
println()

# next, jsave to JSON, process every entry in coupling_info folder.