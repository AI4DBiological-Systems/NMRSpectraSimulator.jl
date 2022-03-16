
# different parameterizations of the spin system FID parameters.
mutable struct SpinSysFIDType1{T}
    κs_λ::Vector{T} # a common multiplier for each spin group. length: number of spin groups.
    κs_β::Vector{Vector{T}} # a vector coefficient for each (spin group). vector length: number of spins in the spin group.
    d::Vector{T} # a common multiplier for each spin group. length: number of spin groups.
end

function SpinSysFIDType1(x::T) where T
    return SpinSysFIDType1{T}(Vector{T}(undef,0), Vector{Vector{T}}(undef, 0), Vector{T}(undef, 0))
end


mutable struct SpinSysFIDType2{T}
    κs_λ::Vector{Vector{T}} # a multiplier for each (spin group, partition element).
    κs_β::Vector{Vector{T}} # a vector coefficient for each (spin group). vector length: number of spins in the spin group.
    d::Vector{Vector{T}} # a multiplier for each (spin group, partition element).
end

function SpinSysFIDType2(x::T) where T
    return SpinSysFIDType2(Vector{Vector{T}}(undef,0), Vector{Vector{T}}(undef, 0), Vector{Vector{T}}(undef, 0))
end

function constructorSSFID(x::SpinSysFIDType1{T}, y...)::SpinSysFIDType1{T} where T
    return SpinSysFIDType1(y...)
end

function constructorSSFID(x::SpinSysFIDType2{T}, y...)::SpinSysFIDType2{T} where T
    return SpinSysFIDType2(y...)
end

###### system-level container.

mutable struct CompoundFIDType{T,SST}

    # resonance components in spin systems.
    qs::Vector{Vector{Function}} # spin group, partition element index.
    qs0::Vector{Vector{Function}} # qs but no phase evaluation. For calibration's least squares.

    αs::Vector{Vector{T}}
    Ωs::Vector{Vector{T}}
    part_inds_compound::Vector{Vector{Vector{Int}}}
    Δc_m_compound::Vector{Vector{Vector{T}}}

    # κs_λ::Vector{T} # multipliers for each spin group.
    # κs_β::Vector{Vector{T}} # coefficients for each (spin group, partition element).
    # d::Vector{T}
    ss_params::SST

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

mutable struct κCompoundFIDType{T,SST}
    κ::Vector{Vector{T}} # spin group, partition element index.
    κ_singlets::Vector{T}
    core::CompoundFIDType{T,SST}
end

function κCompoundFIDType(core::CompoundFIDType{T,SST}) where {T,SST}

    N_spins = length(core.part_inds_compound)
    C = Vector{Vector{T}}(undef, N_spins)
    for i = 1:N_spins
        C[i] = ones(T, length(core.part_inds_compound[i]))
    end

    D = ones(T, length(core.d_singlets))

    return κCompoundFIDType(C, D, core)
end

# creates a reference from As.
function fetchΩS(As::Vector{CompoundFIDType{T,SST}}) where {T,SST}

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
