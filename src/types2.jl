
# different parameterizations of the spin system FID parameters.


struct SpinSysParamsType2{T}
    κs_λ::Vector{Vector{T}} # a multiplier for each (spin group, partition element).
    κs_β::Vector{Vector{T}} # a vector coefficient for each (spin group). vector length: number of spins in the spin group.
    d::Vector{Vector{T}} # a multiplier for each (spin group, partition element).
end

function SpinSysParamsType2(x::T) where T
    return SpinSysParamsType2(Vector{Vector{T}}(undef,0), Vector{Vector{T}}(undef, 0), Vector{Vector{T}}(undef, 0))
end

struct SpinSysParamsType1{T}
    κs_λ::Vector{T} # a common multiplier for each spin group. length: number of spin groups.
    κs_β::Vector{Vector{T}} # a vector coefficient for each (spin group). vector length: number of spins in the spin group.
    d::Vector{T} # a common multiplier for each spin group. length: number of spin groups.
end

function SpinSysParamsType1(x::T) where T
    return SpinSysParamsType1(Vector{T}(undef,0), Vector{Vector{T}}(undef, 0), Vector{T}(undef, 0))
end

function constructorSSFID(x::SpinSysParamsType1{T}, y...)::SpinSysParamsType1{T} where T
    return SpinSysParamsType1(y...)
end

# function constructorSSFID(x::SpinSysFIDType2{T}, y...)::SpinSysFIDType2{T} where T
#     return SpinSysFIDType2(y...)
# end

struct FIDModelType{T,SST} # parameters for surrogate model.

    # non-singlet spin systems.
    qs::Vector{Vector{Function}} # spin group, partition element index.

    # κs_λ::Vector{T} # multipliers for each spin group.
    # κs_β::Vector{Vector{T}} # coefficients for each (spin group, partition element).
    # d::Vector{T}
    ss_params::SST

    # singlets.
    κs_λ_singlets::Vector{T}
    β_singlets::Vector{T}
    d_singlets::Vector{T}

    # misc.
    Δc_max::Vector{T}
    λ0::T
end

mutable struct καFIDModelType{T,SST}
    κs_α::Vector{Vector{T}} # spin group, partition element index.
    κs_α_singlets::Vector{T}
    core::FIDModelType{T,SST}
end

function καFIDModelType(core::FIDModelType{T,SST}) where {T,SST}

    N_spins = length(core.qs)
    κs_α = Vector{Vector{T}}(undef, N_spins)
    for i = 1:N_spins
        κs_α[i] = ones(T, length(core.qs[i]))
    end

    κs_α_singlets = ones(T, length(core.d_singlets))

    return καFIDModelType(κs_α, κs_α_singlets, core)
end

######

struct SHType{T} # output of the SH simulation.

    # resonance components in non-singlet spin systems.
    αs::Vector{Vector{T}}
    Ωs::Vector{Vector{T}}

    Δc_m_compound::Vector{Vector{Vector{T}}}
    part_inds_compound::Vector{Vector{Vector{Int}}}
    Δc_bar::Vector{Vector{Vector{T}}}

    N_spins_sys::Vector{Int}

    # resonance components in singlets spin systems.
    αs_singlets::Vector{T}
    Ωs_singlets::Vector{T}

    # misc input for the SH simulation.
    fs::T
    SW::T
    ν_0ppm::T
end



"""
Remove `Δc_m_compound` and `part_inds_compound` to reduce the memory footprint of Vector{CompoundFIDType} types.
"""
function removeauxinfo!(As::Vector{SHType{T}}) where T

    for n = 1:length(As)
        A = As[n]
        for i = 1:length(A.Δc_m_compound)
            resize!(A.Δc_m_compound[i], 0)
            resize!(A.part_inds_compound[i], 0)
        end
    end

    return nothing
end

function varinfocompound(A)

    property_names = propertynames(A)

    for j = 1:length(property_names)
        display(property_names[j])
        asdf = getfield(A, property_names[j]);

        info_string = varinfo(r"asdf")
        display(info_string)
    end
end
