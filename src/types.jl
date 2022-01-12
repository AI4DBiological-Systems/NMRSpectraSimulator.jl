
mutable struct CompoundFIDType{T}

    # resonance components in spin systems.
    qs::Vector{Vector{Function}} # spin group, partition element index.
    αs::Vector{Vector{T}}
    Ωs::Vector{Vector{T}}
    part_inds_compound::Vector{Vector{Vector{Int}}}
    Δc_m_compound::Vector{Vector{Vector{T}}}

    κs_λ::Vector{T} # multipliers for each spin group.
    κs_β::Vector{Vector{T}} # coefficients for each (spin group, partition element).
    d::Vector{T}

    # resonance components in singlets spin systems.
    αs_singlets::Vector{T}
    Ωs_singlets::Vector{T}

    κs_λ_singlets::Vector{T}
    β_singlets::Vector{T}
    d_singlets::Vector{T}

    # misc info about this compound, and info about the machine settings.
    λ0::T
    fs::T
    SW::T
    Δcs_max::T
    ν_0ppm::T
end

mutable struct κCompoundFIDType{T}
    κ::Vector{Vector{T}} # spin group, partition element index.
    κ_singlets::Vector{T}
    core::CompoundFIDType{T}
end

#const AllCompoundFIDType{T} = Union{CompoundFIDType{T}, κCompoundFIDType{T}}

function κCompoundFIDType(core::CompoundFIDType{T}) where T

    N_spins = length(core.part_inds_compound)
    C = Vector{Vector{T}}(undef, N_spins)
    for i = 1:N_spins
        C[i] = ones(T, length(core.part_inds_compound[i]))
    end

    D = ones(T, length(core.d_singlets))

    return κCompoundFIDType(C, D, core)
end

# creates a reference from As.
function fetchΩS(As::Vector{CompoundFIDType{T}}) where T

    ΩS = Vector{Vector{Vector{T}}}(undef, 0)
    j = 0

    for n = 1:length(As)

        push!(ΩS, As[n].Ωs)
        j += 1

        for i = 1:length(As[n].Ωs_singlets)
            push!(ΩS[j], [ As[n].Ωs_singlets[i]; ])
        end
    end

    return ΩS
end
