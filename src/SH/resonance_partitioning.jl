


#####

function partitionresonances(x_in::Vector{Vector{T}},
    amplitudes_in::Vector{T},
    threshold_amplitude::T;
    radius = 1e-1,
    Minkowski_parameter = 3.5) where T <: Real

    @assert length(x_in) == length(amplitudes_in)

    # take largest amplitude, set as new partition centre.
    amplitudes = copy(amplitudes_in)
    x = copy(x_in)
    inds_buffer = collect(1:length(x))

    out_inds = Vector{Vector{Int}}(undef, 0)

    max_amplitude, max_ind = findmax(amplitudes)
    while !isempty(amplitudes) || max_amplitude < threshold_amplitude

        # set up tree.
        X = Matrix{T}(undef, length(x[1]), length(x))
        for r = 1:length(x)
            X[:,r] = x[r]
        end

        balltree = NearestNeighbors.BallTree(X,
            NearestNeighbors.Minkowski(Minkowski_parameter);
            reorder = false)

        inds = NearestNeighbors.inrange(balltree, x[max_ind], radius, true)

        # book keep.
        push!(out_inds, inds_buffer[inds])

        deleteat!(amplitudes, inds)
        deleteat!(x, inds)
        deleteat!(inds_buffer, inds)

        if isempty(amplitudes)
            return out_inds
        end
        max_amplitude, max_ind = findmax(amplitudes)
    end

    return out_inds
end

"""
αs and Ωs must not contain singlet groups.
"""
function partitionresonances(coherence_state_pairs_sys, ms_sys,
    αs::Vector{Vector{T}}, Ωs::Vector{Vector{T}};
    α_relative_threshold = 0.05) where T

    N_groups = length(αs)

    part_inds_set = Vector{Vector{Vector{Int}}}(undef, N_groups)
    Δc_m_set = Vector{Vector{Vector{T}}}(undef, N_groups)

    # as = Vector{Vector{Vector{T}}}(undef, N_groups)
    # Fs = Vector{Vector{Vector{T}}}(undef, N_groups)
    as = Vector{Vector{T}}(undef, N_groups)
    Fs = Vector{Vector{T}}(undef, N_groups)

    N_spin_sys = length(coherence_state_pairs_sys)
    for i = 1:N_spin_sys

        α_tol = α_relative_threshold*maximum(αs[i])
        inds_amp = findall(xx->(xx>α_tol), αs[i])

        c_states_prune = (coherence_state_pairs_sys[i])[inds_amp]
        αs_i_prune = αs[i][inds_amp]
        Ωs_i_prune = Ωs[i][inds_amp]

        c_m_r = collect( ms_sys[i][r] for (r,s) in c_states_prune )
        c_m_s = collect( ms_sys[i][s] for (r,s) in c_states_prune )
        Δc_m = collect( c_m_r[j] - c_m_s[j] for j = 1:length(c_m_r))

        part_inds = partitionresonances(Δc_m,
            αs_i_prune, α_tol; radius = 1e-1)

        # partition_size = length(part_inds)
        # as[i] = Vector{Vector{T}}(undef, partition_size)
        # Fs[i] = Vector{Vector{T}}(undef, partition_size)
        # for m = 1:partition_size
        #
        #     as[i][m] = αs_i_prune[part_inds[m]]
        #     Fs[i][m] = Ωs_i_prune[part_inds[m]]
        # end
        as[i] = αs_i_prune
        Fs[i] = Ωs_i_prune

        part_inds_set[i] = part_inds
        Δc_m_set[i] = Δc_m
    end

    return as, Fs, part_inds_set, Δc_m_set
end

function fitproxies!(As::Vector{CompoundFIDType{T}};
    κ_λ_lb = 0.5,
    κ_λ_ub = 2.5,
    u_min = Inf,
    u_max = Inf,
    Δr = 1.0,
    Δκ_λ = 0.05) where T

    for n = 1:length(As)
        fitproxy!(As[n];
        κ_λ_lb = κ_λ_lb,
        κ_λ_ub = κ_λ_ub,
        u_min = u_min,
        u_max =u_max,
        Δr = Δr,
        Δκ_λ = Δκ_λ)
    end

    return nothing
end

function fitproxy!(A::CompoundFIDType{T};
    κ_λ_lb = 0.5,
    κ_λ_ub = 2.5,
    u_min = Inf,
    u_max = Inf,
    Δr = 1.0,
    Δκ_λ = 0.05) where T

    hz2ppmfunc = uu->(uu - A.ν_0ppm)*A.SW/A.fs
    ppm2hzfunc = pp->(A.ν_0ppm + pp*A.fs/A.SW)

    # threshold and partition the resonance components.
    if !isfinite(u_min) || !isfinite(u_max)
        Ωs_ppm = hz2ppmfunc.( combinevectors(A.Ωs) ./ (2*π) )
        min_ppm = minimum(Ωs_ppm) - 0.5
        max_ppm = maximum(Ωs_ppm) + 0.5
        u_min = ppm2hzfunc(min_ppm)
        u_max = ppm2hzfunc(max_ppm)
    end

    d_max = ppm2hzfunc(A.Δcs_max)-ppm2hzfunc(0.0)
    A.qs = setupcompoundpartitionitp(d_max,
        A.Δc_m_compound,
        A.part_inds_compound,
        A.αs, A.Ωs,
        A.λ0, u_min, u_max;
        κ_λ_lb = κ_λ_lb,
        κ_λ_ub = κ_λ_ub,
        Δr = Δr,
        Δκ_λ = Δκ_λ)

    return nothing
end

function setupmixtureproxies(target_names::Vector{String},
    base_path_JLD, Δcs_max_mixture::Vector{T},
    hz2ppmfunc, ppm2hzfunc, fs, SW, λ0,
    ν_0ppm::T;
    tol_coherence = 1e-2,
    α_relative_threshold = 0.05) where T <: Real

    N_compounds = length(target_names)
    As = Vector{CompoundFIDType{T}}(undef, N_compounds)

    for n = 1:N_compounds
        As[n] = setupcompoundproxy(target_names[n],
        base_path_JLD, Δcs_max_mixture[n], hz2ppmfunc, ppm2hzfunc, fs, SW, λ0, ν_0ppm;
        tol_coherence = tol_coherence,
        α_relative_threshold = α_relative_threshold)
    end

    return As
end

function setupcompoundproxy(name, base_path, Δcs_max::T, hz2ppmfunc, ppm2hzfunc,
    fs::T, SW::T, λ0::T, ν_0ppm::T;
    tol_coherence = 1e-2,
    α_relative_threshold = 0.05) where T <: Real

    # fetch GISSMO entry.
    records = GISSMOReader.getGISSMOentriesall()
    record_names = collect( records[i].molecule_name for i = 1:length(records) )

    k = findfirst(xx->xx==name, record_names)

    N_targets = length(record_names)
    tmp = collect( GISSMOReader.getGISSMOentry(record_names[i]) for i = 1:N_targets )
    entries = [ records[k].entry ]
    molecule_names = [records[k].molecule_name;]

    record_name = entries[1]
    molecule_name = molecule_names[1]

    # SH.
    p_cs_sys, css_sys, cs_singlets, J_vals_sys, J_IDs_sys, intermediates_sys,
    cs_LUT, cs_singlets_compact, N_spins_singlet, css, p_compound,
    cs_compound = fetchSHparameters(molecule_name, record_name, base_path)

    αs_inp, Ωs_inp, coherence_mat_sys, eig_vals_sys, Q_sys,
    coherence_state_pairs_sys, H_sys, states_sys, ms_sys,
    M_sys = setupcompoundspectrum!(css_sys,
    p_cs_sys, J_vals_sys, J_IDs_sys, intermediates_sys, cs_LUT,
    ppm2hzfunc, cs_singlets_compact, N_spins_singlet, fs, SW;
    tol_coherence = tol_coherence)



    k = findfirst(xx->length(xx)==1, αs_inp)
    αs_spin_sys = copy(αs_inp)
    Ωs_spin_sys = copy(Ωs_inp)
    αs_singlets = Vector{T}(undef, 0)
    Ωs_singlets = Vector{T}(undef, 0)
    if typeof(k) == Int
        αs_spin_sys = αs_spin_sys[1:k-1]
        Ωs_spin_sys = Ωs_spin_sys[1:k-1]

        αs_singlets = collect( αs_inp[l][1] for l = k:length(αs_inp) )
        Ωs_singlets = collect( Ωs_inp[l][1] for l = k:length(Ωs_inp) )
    end

    αs, Ωs, part_inds_compound,
    Δc_m_compound = partitionresonances(coherence_state_pairs_sys,
    ms_sys, αs_spin_sys, Ωs_spin_sys;
    α_relative_threshold = α_relative_threshold)

    # proxy placeholder.
    qs = Vector{Vector{Function}}(undef, 0)

    # spectra eval func for spin systems.
    N_spins_compound = collect( length(css_sys[i]) for i = 1:length(css_sys))
    #κs_λ = collect( convertcompactdomain(rand(), 0.0, 1.0, κ_λ_lb, κ_λ_ub) for i = 1:length(Ωs_inp))
    #κs_β = collect( randn(N_spins_compound[i]) for i = 1:length(N_spins_compound))
    κs_λ = ones(T, length(Ωs_inp))
    κs_β = collect( zeros(T, N_spins_compound[i]) for i = 1:length(N_spins_compound))
    #d = rand(length(αs))
    d = zeros(T, length(αs))

    N_singlets = length(αs_singlets)
    κs_λ_singlets = ones(T, N_singlets)
    κs_β_singlets = zeros(T, N_singlets)
    d_singlets = zeros(T, N_singlets)

    # f = uu->evalcLcompoundviapartitions(uu, d,
    # αs, Ωs, κs_λ, κs_β, λ0, Δc_m_compound, part_inds_compound)

    out = CompoundFIDType(qs, αs, Ωs, part_inds_compound, Δc_m_compound,
    κs_λ, κs_β, d,
    αs_singlets, Ωs_singlets, κs_λ_singlets, κs_β_singlets, d_singlets,
    λ0, fs, SW, Δcs_max, ν_0ppm)

    return out
end


function fetchsubset(As::Vector{CompoundFIDType{T}}, indices::Vector{Tuple{Int,Int}}) where T <: Real

    # get
    inds = collect( findall(xx->xx[1]==n, indices) for n = 1:length(As) )
    filter!(xx->!isempty(xx), inds)

    Bs = Vector{CompoundFIDType{T}}(undef, length(inds))

    for j = 1:length(inds)
        ks = inds[j]
        n, _ = indices[ks[1]]
        A = As[n]

        is = collect( indices[ks[l]][2] for l = 1:length(ks) )
        is_sys = filter(xx->xx<=length(A.qs), is)

        N_sys = length(A.qs)
        is_singlets = filter(xx->xx>length(A.qs), is) .- N_sys


        Bs[j] = CompoundFIDType(A.qs[is_sys], A.αs[is_sys], A.Ωs[is_sys],
        A.part_inds_compound[is_sys], A.Δc_m_compound[is_sys],
        A.κs_λ[is_sys], A.κs_β[is_sys], A.d[is_sys],
        A.αs_singlets[is_singlets],
        A.Ωs_singlets[is_singlets],
        A.κs_λ_singlets[is_singlets],
        A.β_singlets[is_singlets],
        A.d_singlets[is_singlets],
        A.λ0, A.fs, A.SW, A.Δcs_max, A.ν_0ppm)
    end

    return Bs
end

# I am here.
function getNβ(A::CompoundFIDType{T}) where T

    counter_sys = 0
    for i = 1:length(A.κs_β)
        counter_sys += length(A.κs_β[i])
    end

    N_β = counter_sys + length(A.β_singlets)

    return N_β
end
