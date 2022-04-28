
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
    gs::Vector{Vector{Function}} # qs but no phase evaluation. For calibration's least squares.

    αs::Vector{Vector{T}}
    Ωs::Vector{Vector{T}}
    part_inds_compound::Vector{Vector{Vector{Int}}}
    Δc_m_compound::Vector{Vector{Vector{T}}}

    Δc_bar::Vector{Vector{Vector{T}}}

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
    ν_0ppm::T
end

mutable struct κCompoundFIDType{T,SST}
    κs_α::Vector{Vector{T}} # spin group, partition element index.
    κs_α_singlets::Vector{T}
    core::CompoundFIDType{T,SST}
end


"""
Remove gs and Δc_m_compound to reduce the memory footprint of Vector{CompoundFIDType} types.
"""
function removeauxinfo!(As::Vector{CompoundFIDType{T,SST}}) where {T,SST}

    for n = 1:length(As)
        A = As[n]
        for i = 1:length(A.Δc_m_compound)
            resize!(A.Δc_m_compound[i], 0)
        end

        for i = 1:length(A.gs)
            resize!(A.gs[i], 0)
        end
    end

    return nothing
end

function varinfocompound(A::CompoundFIDType{T,SST}) where {T,SST}

    property_names = propertynames(A)
    for j = 1:length(property_names)
        display(property_names[j])
        asdf = getfield(A, property_names[j]);

        info_string = varinfo(r"asdf")
        display(info_string)
    end
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

mutable struct zCompoundFIDType{T,SST}
    zs::Vector{Vector{Complex{T}}} # spin group, partition element index.
    zs_singlets::Vector{Complex{T}}
    core::CompoundFIDType{T,SST}
end

function zCompoundFIDType(core::CompoundFIDType{T,SST}) where {T,SST}

    N_spins = length(core.part_inds_compound)
    C = Vector{Vector{Complex{T}}}(undef, N_spins)
    for i = 1:N_spins
        C[i] = ones(Complex{T}, length(core.part_inds_compound[i]))
    end

    D = ones(Complex{T}, length(core.d_singlets))

    return zCompoundFIDType(C, D, core)
end
