
import JSON
import JSON3

# assigns a default value for some entries for every key (compound name) in dict_compound_to_filename.

save_folder = "/home/roy/Documents/repo/NMRData/input/SH_configs"
save_name = "select_compounds_SH_config.json"
save_path = joinpath(save_folder, save_name)

dict_compound_to_filename = JSON.parsefile("/home/roy/Documents/repo/NMRData/input/compound_mapping/select_compounds.json")

function createdefaultSHconfig(guiding_dict)
    db_dict = Dict()
    for (key, value) in guiding_dict

        db_dict[key] = Dict("coherence tolerance" => 1e-2,
            "relative amplitude threshold" => 0.05,
            "maximum Δc deviation" => 1e-1)
    end

    return db_dict
end

dict_out = createdefaultSHconfig(dict_compound_to_filename)

stringdata = JSON.json(dict_out)

open(save_path, "w") do f
    JSON3.pretty(f, stringdata)
    println(f)
end
