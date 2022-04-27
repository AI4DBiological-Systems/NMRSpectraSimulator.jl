function addcompoundtoconfig!(out_string::Vector{String},
    molecule_name::String,
    css_sys::Vector{Vector{T}},
    cs_singlets_compact::Vector{T},
    attribute_text::Vector{String},
    default_values) where T <: Real

    # save heading to template.
    push!(out_string, molecule_name)

    # save the default cs configuration info.
    N_spin_sys = length(css_sys)
    for i = 1:N_spin_sys
        push!(out_string, "Spin group $(i)")

        for k = 1:length(attribute_text)
            line_string = "$(attribute_text[k]) $(default_values[k])"
            push!(out_string, line_string)
        end


        push!(out_string, "end group")
    end

    for i = 1:length(cs_singlets_compact)
        push!(out_string, "Singlet $(i)")

        for k = 1:length(attribute_text)
            line_string = "$(attribute_text[k]) $(default_values[k])"
            push!(out_string, line_string)
        end

        push!(out_string, "end singlet")
    end


    push!(out_string, "end compound")
    push!(out_string, "")

    return nothing
end



function createdefaultconfig(molecule_names,
    css_sys, cs_singlets_compact,
    base_dir::String,
    save_dir::String;
    δ_lb::T = 0.1,
    δ_ub::T = 0.1,
    unique_cs_tol::T = 1e-6,
    delay_time = 0.2,
    default_cs_delta::T = 0.1) where T <: Real

    N = length(entries)

    out_string = Vector{String}(undef, 0)

    for i = 1:N

        addentrydefaultconfigtobuffer!(out_string,
            entries[i].molecule_name,
            css_sys,
            cs_singlets_compact;
            default_cs_delta = default_cs_delta)
    end

    # special cases: add D2O.
    entry = "0"
    molecule_name = "D2O"
    singlet_ppm = 4.8
    css_sys, cs_singlets_compact = specialsingletentry(entry,
        molecule_name,
        base_dir,
        save_dir,
        singlet_ppm;
            δ_lb = δ_lb,
            δ_ub = δ_ub,
            unique_cs_tol = unique_cs_tol)

    # add to buffer.
    addentrydefaultconfigtobuffer!(out_string,
            molecule_name,
            css_sys,
            cs_singlets_compact;
            default_cs_delta = default_cs_delta)

    # write buffer to file.
    file_path = joinpath(save_dir, "default_cs_config.txt")
    DelimitedFiles.writedlm(file_path, out_string)

    return nothing
end
