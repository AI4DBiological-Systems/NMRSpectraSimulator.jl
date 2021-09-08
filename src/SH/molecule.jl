# for multiple molecules.
function setupfit1molecules(    base_path::String,
    entries::Vector{String},
    molecule_names::Vector{String},
    ppm2hzfunc,
    fs::T,
    SW::T;
    tol_coherence = 1e-2, # 1e-7 #for -1 coherence.
    prune_a_tol_factor = 0.01,
    prune_zero_tol = 1e-8,
    fuse_radius = ppm2hzfunc(0.001) - ppm2hzfunc(0.0),
    Minkowski_parameter = 3.5,
    partition_radius = ppm2hzfunc(0.1)-ppm2hzfunc(0.0),
    partition_Minkowski_parameter = 3.5,
    update_α_flag::Bool = false) where T <: Real

    N_compounds = length(molecule_names)
    N = N_compounds

    updateαΩFSfunc_set = Vector{Function}(undef, N)

    # approximation output.
    αs_set = Vector{Vector{Vector{T}}}(undef, N)
    Ωs_set = Vector{Vector{Vector{T}}}(undef, N)

    # cs.
    p0_set = Vector{Vector{T}}(undef, N)
    N_p_cs_sys_set = Vector{Vector{Int}}(undef, N)
    p0_group_set::Vector{Vector{Vector{T}}} = Vector{Vector{Vector{T}}}(undef, N)

    # compounds.
    for n = 1:N_compounds
    αs_set[n], Ωs_set[n],
    updateαΩFSfunc_set[n],
    p0_set[n],
    N_p_cs_sys_set[n], p0_group_set[n] = setupfit1molecule(base_path,
            entries[n], molecule_names[n],
            ppm2hzfunc,
            fs,
            SW;
            tol_coherence = tol_coherence,
            prune_a_tol_factor = prune_a_tol_factor,
            prune_zero_tol = prune_zero_tol,
            fuse_radius = fuse_radius,
            Minkowski_parameter = Minkowski_parameter,
            partition_radius = partition_radius,
            partition_Minkowski_parameter = partition_Minkowski_parameter,
            update_α_flag = update_α_flag)
    end

    return αs_set, Ωs_set,
    updateαΩFSfunc_set, p0_set,
    N_p_cs_sys_set, p0_group_set
end

function setupfit1molecule(  base_path,
    record_name,
    molecule_name,
    ppm2hzfunc,
    fs,
    SW;
    tol_coherence = 1e-2, # 1e-7 #for -1 coherence.
    prune_a_tol_factor = 0.05,
    prune_zero_tol = 1e-8,
    fuse_radius = ppm2hzfunc(0.001) - ppm2hzfunc(0.0),
    Minkowski_parameter = 3.5,
    partition_radius = ppm2hzfunc(0.1)-ppm2hzfunc(0.0),
    partition_Minkowski_parameter = 3.5,
    update_α_flag::Bool = false) where T <: Real

    # load molecule information.
    file_name = "$(molecule_name)_$(record_name)"

    # tmp = filter(x->occursin(molecule_name,x), readdir(base_path))
    # file_name = tmp[1][1:end-4] # take first match, and get rid of the extension ".jld".

    css_sys, J_IDs_sys, J_vals_sys, cs_singlets,
    J_lb_sys, J_ub_sys, cs_lb_sys, cs_ub_sys,
    cs_LUT, reference_concentration, solute_concentration, _ = loadGISSMOmolecule(base_path, file_name)

    #N_groups = length(cs_group)
    p_cs_sys = cstopcs(css_sys, cs_LUT)

    #cs_singlets_compact = collect( cs_singlets[m][1] for m = 1:length(cs_singlets) )
    #N_spins_singlet = collect( length(cs_singlets[m]) for m = 1:length(cs_singlets) )

    N_spins_singlet = Vector{Int}(undef, length(cs_singlets))
    cs_singlets_compact = Vector{Float64}(undef, length(cs_singlets))
    for m = 1:length(cs_singlets)
        cs_singlets_compact[m] = cs_singlets[m][1]
        N_spins_singlet[m] = length(cs_singlets[m])
    end

    N_sys = length(css_sys)
    @assert length(p_cs_sys) == N_sys

    N_spins_sys = collect( length(css_sys[m]) for m = 1:length(css_sys) )
    intermediates_sys = prepcouplingalgorithm(N_spins_sys)

    # prepare initial guess and bounds.
    #tmp = collect( [cs_group[i]; (J_lb_group[i] + J_ub_group[i]) ./ 2] for i = 1:N_groups )
    css = combinevectors(p_cs_sys)
    push!(css, cs_singlets_compact...)

    N_p_cs_sys = collect( length(p_cs_sys[i]) for i = 1:N_sys )
    p0 = copy(css)

    # could be used for dim-reduced p.
    p0_group = Vector{Vector{Float64}}(undef, N_sys + length(cs_singlets_compact))
    for i = 1:N_sys
        #p0_group[i] = p_cs_sys[i]
        p0_group[i] = copy(p_cs_sys[i])
    end
    for i = 1:length(cs_singlets_compact)
        p0_group[N_sys+i] = [ cs_singlets_compact[i] ]
    end

    #### set up.
    αs, eig_vals, updateΩFSfuncs, Ωs,
    states, state_labels = setupapproximateΩmolecule!( css_sys,
    p_cs_sys,
    J_vals_sys,
    J_IDs_sys,
    intermediates_sys,
    cs_LUT,
    ppm2hzfunc,
    cs_singlets_compact,
    N_spins_singlet,
    N_spins_sys,
    fs,
    SW;
    tol_coherence = tol_coherence,
    prune_a_tol_factor = prune_a_tol_factor,
    fuse_radius = fuse_radius,
    Minkowski_parameter = Minkowski_parameter,
    update_α_flag = update_α_flag)

    ##### strong coupling.

    ## strong coupling.
    cs_len_sys = getcslengthfromLUT(cs_LUT)

    # Full solve. Takes Δp_cs as input.
    # separate buffer for css_sys, since we're storing Δcs.
    N_singlets = length(cs_singlets_compact)
    cs_singlets_compact0 = copy(cs_singlets_compact)
    updateαΩFSfunc = (pp,ii)->evalfullsolveproxy!(deepcopy(css_sys),
                cs_singlets_compact0,
                #N_singlets,
                Ωs,
                pp,
                ii,
                cs_LUT,
                updateΩFSfuncs,
                fs, SW, ppm2hzfunc)


    return αs, Ωs, updateαΩFSfunc,
    p0, N_p_cs_sys, p0_group
end

function setupfit1experiment(experiment_full_path::String,
base_path_JLD::String,
entries_input::Vector{String},
molecule_names_input::Vector{String},
ref_entry::String,
ref_name::String,
solvent_entry::String,
solvent_name::String;
tol_coherence = 1e-2, # 1e-7 #for -1 coherence.
prune_a_tol_factor = 0.01,
fuse_radius_δppm = 0.001,
partition_radius_δppm = 0.1,
prune_zero_tol = 1e-8,
Minkowski_parameter = 3.5,
partition_Minkowski_parameter = 3.5,
add_reference_and_solvent_flag = true,
solvent_ppm_guess = 4.8,
solvent_window_ppm = 0.1,
update_α_flag::Bool = false)


    # load the experiment.
    s_t, S, hz2ppmfunc, ppm2hzfunc, ν_0ppm, fs, SW,
    α_0ppm, β_0ppm, λ_0ppm, Ω_0ppm,
       results_0ppm, dic, α_solvent, β_solvent, λ_solvent, Ω_solvent,
        results_solvent = loadspectrum(experiment_full_path;
                                    solvent_ppm = solvent_ppm_guess,
                                    solvent_window_ppm = solvent_window_ppm)

    # set up SC and proxy parameters.
    fuse_radius = ppm2hzfunc(fuse_radius_δppm) - ppm2hzfunc(0.0)
    partition_radius = ppm2hzfunc(partition_radius_δppm)-ppm2hzfunc(0.0)

    αs_set, Ωs_set, updateαΩFSfunc_set, p0_set,
    N_p_cs_sys_set, p0_group_set = setupfit1molecules(base_path_JLD,
        entries_input,
        molecule_names_input,
        ppm2hzfunc,
        fs,
        SW;
        tol_coherence = tol_coherence,
        prune_a_tol_factor = prune_a_tol_factor,
        prune_zero_tol = prune_zero_tol,
        fuse_radius = fuse_radius,
        Minkowski_parameter = Minkowski_parameter,
        partition_radius = partition_radius,
        partition_Minkowski_parameter = partition_Minkowski_parameter,
        update_α_flag = update_α_flag)

    return s_t, S, hz2ppmfunc, ppm2hzfunc, ν_0ppm, fs, SW,
        α_0ppm, λ_0ppm, Ω_0ppm, β_0ppm, dic,
        αs_set, Ωs_set,
        updateαΩFSfunc_set, p0_set,
        N_p_cs_sys_set, p0_group_set,
        α_solvent, β_solvent, λ_solvent, Ω_solvent
end
