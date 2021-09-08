

mutable struct ProfilerConfigType{T}
    tol_coherence::T
    prune_zero_tol::T
    partition_radius_δppm::T
    Minkowski_parameter::T
    partition_Minkowski_parameter::T
    prune_a_tol_factor::T
    fuse_radius_δppm::T
    add_reference_and_solvent_flag::Bool
    exclude_0ppm_and_solvent_cs_flag::Bool
    update_α_flag::Bool

    # for choosing test positions.
    N_samples_r::Int
    N_samples_i::Int
    N_samples_abs::Int
    ln_transform_lower::T
    ln_transform_upper_r::T
    ln_transform_upper_i::T
    ln_transform_upper_abs::T
    radius_ppm::T # distance away from ΩS, for pruning fit positions.
    zero_ppm_offset_ppm::T

    solvent_ppm_guess::T
    solvent_window_ppm::T

    # not used in calibration.
    solvent_entry::String
    solvent_name::String
    ref_entry::String
    ref_name::String
end

# add_reference_and_solvent_flag = false for calibration.
function defaultconfig( add_reference_and_solvent_flag::Bool,
    update_α_flag::Bool)

    #L_min = 8
    tol_coherence = 1e-2 # 1e-7 #for -1 coherence.

    exclude_0ppm_and_solvent_cs_flag = false # consider removing this option.

    prune_zero_tol = 1e-8
    partition_radius_δppm = 0.1
    Minkowski_parameter = 3.5
    partition_Minkowski_parameter = 3.5

    prune_a_tol_factor = 0.01
    fuse_radius_δppm = 0.001

    #prune_a_tol_factor = 0.001
    #fuse_radius_δppm = 0.0005
    #fuse_radius_δppm = 0.000

    N_samples_r = 2000
    N_samples_i = 2000
    N_samples_abs = 2000
    # N_samples_r = 4000
    # N_samples_i = 4000
    # N_samples_abs = 4000
    ln_transform_lower = 1.0
    ln_transform_upper_r = 10.0
    ln_transform_upper_i = 10.0
    ln_transform_upper_abs = 30.0

    radius_ppm = 0.3
    zero_ppm_offset_ppm = -0.3

    solvent_ppm_guess = 4.8
    solvent_window_ppm = 0.1

    solvent_entry = "0"
    solvent_name = "D2O"

    tmp = getGISSMOentry("DSS")
    ref_entry = tmp.entry
    ref_name = tmp.molecule_name



    config = ProfilerConfigType(tol_coherence,
        prune_zero_tol,
        partition_radius_δppm,
        Minkowski_parameter,
        partition_Minkowski_parameter,
        prune_a_tol_factor,
        fuse_radius_δppm,
        add_reference_and_solvent_flag,
        exclude_0ppm_and_solvent_cs_flag,
        update_α_flag,

        # for choosing test positions.
        N_samples_r,
        N_samples_i,
        N_samples_abs,
        ln_transform_lower,
        ln_transform_upper_r,
        ln_transform_upper_i,
        ln_transform_upper_abs,
        radius_ppm,
        zero_ppm_offset_ppm,

        solvent_ppm_guess,
        solvent_window_ppm,

        # not used in calibration.
        solvent_entry,
        solvent_name,
        ref_entry,
        ref_name)
    return config
end

# calibrate to a particular GISSMO simulation.
function runcalibration(experiment_full_path::String,
    base_path_JLD::String,
    entries::Vector{String},
    molecule_names::Vector{String},
    config::ProfilerConfigType{T}) where T

    if config.add_reference_and_solvent_flag
    # add reference compound.
    push!(entries, config.ref_entry)
    push!(molecule_names, config.ref_name)

    # add solvent compound.
    push!(entries, config.solvent_entry)
    push!(molecule_names, config.solvent_name)
    end

    # set up p_cs proxies.
    s_t, S, hz2ppmfunc, ppm2hzfunc, ν_0ppm, fs, SW,
    α_0ppm, λ_0ppm, Ω_0ppm, β_0ppm, data_dic,
    as_set, Fs_set, updateαΩFSfunc_set, p0_set,
    N_p_cs_sys_set, p0_group_set,
    α_solvent, β_solvent, λ_solvent,
    Ω_solvent = setupfit1experiment( experiment_full_path,
        base_path_JLD,
        entries,
        molecule_names,
        config.ref_entry,
        config.ref_name,
        config.solvent_entry,
        config.solvent_name;
        tol_coherence = config.tol_coherence,
        prune_a_tol_factor = config.prune_a_tol_factor,
        fuse_radius_δppm = config.fuse_radius_δppm,
        partition_radius_δppm = config.partition_radius_δppm,
        prune_zero_tol = config.prune_zero_tol,
        Minkowski_parameter = config.Minkowski_parameter,
        partition_Minkowski_parameter = config.partition_Minkowski_parameter,
        add_reference_and_solvent_flag = config.add_reference_and_solvent_flag,
        solvent_ppm_guess = config.solvent_ppm_guess,
        solvent_window_ppm = config.solvent_window_ppm,
        update_α_flag = config.update_α_flag)

    N_compounds = length(p0_set)
    #cs_set = collect( sort(p0_set[n]) for n = 1:N_compounds)
    N_cs_vars = sum( length(p0_set[n]) for n = 1:N_compounds)
    N_spin_groups = collect( length(p0_set[n]) for n = 1:N_compounds )

    # force combine of fullsolve.
    #println("updateαΩFSfunc_set = ", updateαΩFSfunc_set)

    p_cs0 = combinevectors(p0_set)
    Δp_cs = zeros(T, length(p_cs0))
    cs_st_ind = 1
    cs_fin_ind = updatemixtureαΩ( updateαΩFSfunc_set, Δp_cs, cs_st_ind)
    ΩS_FS = deepcopy(Fs_set)
    αS_FS = deepcopy(as_set)

    #return αS_FS, ΩS_FS, αS_Rayleigh, ΩS_Rayleigh, as_set, Fs_set

    # # sanity-check.
    # @assert norm(combinevectors(combinevectors(ΩS_Rayleigh)) -
    # combinevectors(combinevectors(ΩS_FS)) ) < 1e-3
    #
    # if !config.update_α_flag
    # @assert norm(combinevectors(combinevectors(αS_Rayleigh)) -
    # combinevectors(combinevectors(αS_FS)) ) < 1e-5
    # end

    #return αS_FS, ΩS_FS, αS_Rayleigh, ΩS_Rayleigh, as_set, Fs_set

    return s_t, S, hz2ppmfunc, ppm2hzfunc, ν_0ppm, fs, SW,
    α_0ppm, λ_0ppm, Ω_0ppm, β_0ppm, data_dic,
    as_set, Fs_set, updateαΩFSfunc_set,
    p0_set, N_p_cs_sys_set, p0_group_set, as_set, Fs_set,
    α_solvent, β_solvent, λ_solvent, Ω_solvent
end


### for updating mixtures.
function updatemixtureαΩ(   updateαΩfunc_set,
                            p::Vector{T},
                            st_ind::Int) where T <: Real

    N_compounds = length(updateαΩfunc_set)

    fin_ind::Int = 0
    for n = 1:N_compounds
        # for i = 1:length(updateαΩfunc_set[n])
        #     st_ind = updateαΩfunc_set[n][i](p, st_ind)
        # end
        fin_ind = updateαΩfunc_set[n](p, st_ind)

        st_ind = fin_ind + 1
    end

    return fin_ind
end
