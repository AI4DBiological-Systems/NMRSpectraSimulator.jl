



function evalcompound(u_rad, A::SHType{T}, B::FIDModelType{T,SST})::Complex{T} where {T <: Real, SST}

    #u_rad = 2*π*u

    out_sys = evalcLcompoundviapartitions(u_rad, A.αs, A.Ωs, B.ss_params,
    B.λ0, A.Δc_bar, A.part_inds_compound)

    out_singlets = evalsinglets(u_rad, B.d_singlets, A.αs_singlets, A.Ωs_singlets,
    B.β_singlets, B.λ0, B.κs_λ_singlets)

    return out_sys + out_singlets
end

function evalmixture(u_rad, As::Vector{SHType{T}}, Bs::Vector{FIDModelType{T,SST}};
    w::Vector{T} = ones(T, length(As)))::Complex{T} where {T <: Real, SST}

    #u_rad = 2*π*u

    out = zero(Complex{T})
    for n = 1:length(As)

        out += w[n]*evalcompound(u_rad, As[n], Bs[n])
    end

    return out
end



"""
Evaluates an element from the partition of a spin group.
r := 2*π*u-d, d is chem shift of this element in radians.
Does not evaluate the complex phase β.
"""
function evalcLpartitionelement(r,
    α::Vector{T}, Ω::Vector{T}, λ::T)::Complex{T} where T <: Real

    out = sum( α[l]/(λ+im*(r-Ω[l])) for l = 1:length(α) )

    return out
end

function evalcLcompoundviapartitions(u_rad,
    αs::Vector{Vector{T}}, Ωs::Vector{Vector{T}},
    x::SpinSysParamsType1{T}, λ0::T,
    c, part_inds_compound)::Complex{T} where T <: Real

    #u_rad = 2*π*u

    out = zero(Complex{T})
    for i = 1:length(αs)
        r = u_rad - x.d[i]

        for k = 1:length(part_inds_compound[i])
            inds = part_inds_compound[i][k]

            out += evalcLpartitionelement(r, αs[i][inds],
                Ωs[i][inds], x.κs_λ[i]*λ0)*exp(im*dot(x.κs_β[i], c[i][k]))
        end
    end

    return out
end

function evalcLcompoundviapartitions(u_rad,
    αs::Vector{Vector{T}}, Ωs::Vector{Vector{T}},
    x::SpinSysParamsType2{T}, λ0::T,
    c, part_inds_compound)::Complex{T} where T <: Real

    #u_rad = 2*π*u

    out = zero(Complex{T})
    for i = 1:length(αs)

        for k = 1:length(part_inds_compound[i])
            r = u_rad - x.d[i][k]
            inds = part_inds_compound[i][k]

            out += evalcLpartitionelement(r, αs[i][inds],
                Ωs[i][inds], x.κs_λ[i][k]*λ0)*exp(im*dot(x.κs_β[i], c[i][k]))
        end
    end

    return out
end
