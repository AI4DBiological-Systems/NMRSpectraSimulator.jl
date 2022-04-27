using LinearAlgebra


import GISSMOReader

include("../src/NMRSpectraSimulator.jl")
import .NMRSpectraSimulator

import Kronecker
import Graphs
import JSON, JSON3

include("./helpers/operators.jl")


H_params_path = "/home/roy/Documents/repo/NMRData/input/coupling_info"
dict_compound_to_filename = JSON.parsefile("/home/roy/Documents/repo/NMRData/input/compound_mapping/select_compounds.json")

name = "L-Histidine"

unique_cs_tol = 1e-6
zero_tol_sigdigits = 6

# load.
load_path = joinpath(H_params_path, dict_compound_to_filename[name]["file name"])
H_IDs, H_css, J_IDs, J_vals = NMRSpectraSimulator.loadcouplinginfojson(load_path)


# process.
#H_inds = collect(1:length(H_IDs))
#J_inds = convertJIDstoJinds(J_IDs, H_IDs)
#g = constructspinsystemsgraphforcompound(length(H_IDs), J_inds)

J_inds_sys2, J_IDs_sys2, J_vals_sys2, H_inds_sys2,
    cs_sys2, H_inds_singlets2, cs_singlets2, H_inds, J_inds,
    g = NMRSpectraSimulator.setupcsJ(H_IDs, H_css, J_IDs, J_vals)

dict_H_inds_to_css = Dict(H_inds .=> H_css)
dict_H_ID_to_ind = Dict(H_IDs .=> collect(1:length(H_IDs)))
dict_ind_to_H_ID = Dict( collect(1:length(H_IDs)) .=> H_IDs)

dict_J_ID_to_val = Dict(J_IDs .=> J_vals)
dict_J_ind_to_val = Dict(J_inds .=> J_vals)

dict_J_ID_to_ind = Dict(J_IDs .=> J_inds)
dict_J_ind_to_ID = Dict(J_inds .=> J_IDs)


systems_g = Graphs.connected_components(g) # node inds for spin systems.
C = Graphs.maximal_cliques(g)

T = Float64
i = 1






### legacy




J_IDs_sys, J_vals_sys, H_IDs_sys = GISSMOReader.partitionJcouplings(J_IDs, J_vals)

css_sys = GISSMOReader.getcssys(H_IDs_sys, H_IDs, H_css)

𝐽_IDs_sys = GISSMOReader.getJsys(J_IDs_sys, J_vals_sys, H_IDs_sys)

H_singlets, cs_singlets,
    cs_singlets_compact = GISSMOReader.createsingletsysems(H_IDs, H_IDs_sys,
                                H_css;
                                zero_tol = unique_cs_tol)
#
### remove spin groups that only have one unique chem shift, but have J-coupling between the atoms.
# Each of such groups are are really one singlet group.
GISSMOReader.removeisolatedcs!(css_sys,
                cs_singlets,
                cs_singlets_compact,
                H_singlets,
                H_IDs_sys,
                J_IDs_sys,
                J_vals_sys,
                𝐽_IDs_sys;
                zero_tol_sigdigits = zero_tol_sigdigits)
#
println("Name: ", name)
println("cs_singlets: ", cs_singlets)
println("cs_singlets_compact: ", cs_singlets_compact)
println("H_singlets: ", H_singlets)

println("css_sys: ", css_sys)
println("J_vals_sys: ", J_vals_sys)
println("𝐽_IDs_sys: ", 𝐽_IDs_sys)
println()
#### end code for moving spin group to singlet group.

#
cs_LUT, p_cs_sys = GISSMOReader.constructLUTcss(css_sys; zero_tol = unique_cs_tol)


cs_len_sys = GISSMOReader.getcslengthfromLUT(cs_LUT)
