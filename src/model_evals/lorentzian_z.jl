# with z-proxy. refactor/remove this later.
function evalitpproxymixture(u, Es::Vector{zCompoundFIDType{T,SST}};
    w::Vector{T} = ones(T, length(Es)))::Complex{T} where {T <: Real, SST}

    u_rad = 2*π*u

    out = zero(Complex{T})

    for n = 1:length(Es)
        out += w[n]*evalitpproxycompound(u, Es[n])
    end

    return out
end

# with z-proxy.
function evalitpproxycompound(u, A::zCompoundFIDType{T,SST})::Complex{T} where {T <: Real, SST}

    u_rad = 2*π*u

    out_sys = evalzitpproxysys(A.zs, A.core.gs, u, A.core.ss_params)

    out_singlets = evalzsinglets(u, A.core.d_singlets,
    A.core.αs_singlets, A.core.Ωs_singlets, A.core.λ0,
    A.core.κs_λ_singlets, A.zs_singlets)

    return out_sys + out_singlets
end

function evalzsinglets(u::T, d::Vector{T},
    αs_singlets::Vector{T}, Ωs_singlets, λ0::T, λ_multipliers::Vector{T};
    zs_singlets::Vector{Complex{T}}) where T <: Real

    u_rad = 2*π*u

    out = zero(Complex{T})
    for i = 1:length(αs_singlets)
        τ = u_rad - d[i]

        λ = λ0*λ_multipliers[i]
        Ω = Ωs_singlets[i]
        out += zs_singlets[i]*αs_singlets[i]/(λ+im*(τ-Ω))
    end
    return out
end

function evalzitpproxysys(zs::Vector{Vector{Complex{T}}}, gs::Vector{Vector{Function}},
    u::T, x::SpinSysFIDType1{T})::Complex{T} where T

    d = x.d
    κs_λ = x.κs_λ

    @assert length(d) == length(gs)

    out = zero(Complex{T})

    u_rad = 2*π*u
    for i = 1:length(gs)
        r = u_rad - d[i]

        for k = 1:length(gs[i])
            out += zs[i][k]*gs[i][k](r, κs_λ[i])
        end
    end

    return out
end

function evalκitpproxysys(zs::Vector{Vector{Complex{T}}}, gs::Vector{Vector{Function}},
    u::T, x::SpinSysFIDType2{T})::Complex{T} where T

    d = x.d
    κs_λ = x.κs_λ

    @assert length(d) == length(gs)

    out = zero(Complex{T})

    u_rad = 2*π*u
    for i = 1:length(gs)

        for k = 1:length(gs[i])
            r = u_rad - d[i][k]

            out += zs[i][k]*gs[i][k](r, κs_λ[i][k])
        end
    end

    return out
end
