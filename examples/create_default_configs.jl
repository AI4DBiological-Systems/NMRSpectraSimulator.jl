# This script create default configuration files that contain the tuning parameters
# for some of the algorithms this package. We assign a default value for each
# compound name key entry in the user-specified
# "compound to parameters filename" JSON file.

# One should manually edit the generated config files according to
# their preferred tuning parameters for the different compounds.

import JSON
import JSON3

name_map_filename = "select_compounds"
load_path = joinpath("/home/roy/Documents/repo/NMRData/input/compound_mapping", "$(name_map_filename).json")
dict_compound_to_filename = JSON.parsefile(load_path)
compound_names = collect( key for (key,value) in dict_compound_to_filename )


##### config file for SH simulation.
function createdefaultSHconfig(compound_names::Vector{String};
    tol_coherence::Float64 = 1e-2,
    α_relative_threshold::Float64 = 0.05,
    Δc_partition_radius::Float64 = 0.3)

    db_dict = Dict()
    for name in compound_names

        db_dict[name] = Dict("coherence tolerance" => tol_coherence,
            "relative amplitude threshold" => α_relative_threshold,
            "maximum Δc deviation" => Δc_partition_radius)
    end

    return db_dict
end

save_folder = "/home/roy/Documents/repo/NMRData/input/SH_configs"
save_name = "$(name_map_filename)_SH_configs_default.json"
save_path = joinpath(save_folder, save_name)

# write the default SH configurations to json.
dict_out = createdefaultSHconfig(compound_names)

stringdata = JSON.json(dict_out)

open(save_path, "w") do f
    JSON3.pretty(f, stringdata)
    println(f)
end


##### config file for SH parameters.
function createdefaultsurrogateconfig(compound_names::Vector{String};
    Δr = 1.0,
    Δκ_λ = 0.05,
    Δc_max::Vector{Float64} = Vector{Float64}(undef, 0),
    κ_λ_lb = 0.5,
    κ_λ_ub = 2.5)

    db_dict = Dict()
    for name in compound_names

        db_dict[name] = Dict("radians surrogate sampling step size" => Δr,
            "κ_λ surrogate sampling step size" => Δκ_λ,
            "Δc_max" => Δc_max,
            "κ_λ lower bound" => κ_λ_lb,
            "κ_λ upper bound" => κ_λ_ub)
    end

    return db_dict
end

save_folder = "/home/roy/Documents/repo/NMRData/input/surrogate_configs"
save_name = "$(name_map_filename)_surrogate_configs_default.json"
save_path = joinpath(save_folder, save_name)

# write the default SH configurations to json.
dict_out = createdefaultsurrogateconfig(compound_names)

stringdata = JSON.json(dict_out)

open(save_path, "w") do f
    JSON3.pretty(f, stringdata)
    println(f)
end
