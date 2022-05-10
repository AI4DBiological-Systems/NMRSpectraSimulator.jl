################### spin system.

function evalitpproxysys(qs::Vector{Vector{Function}},
    u_rad::T, x::SpinSysParamsType1{T})::Complex{T} where T

    d = x.d
    κs_λ = x.κs_λ
    #κs_β = x.κs_β

    @assert length(d) == length(qs)

    out = zero(Complex{T})

    for i = 1:length(qs)
        r = u_rad - d[i]

        for k = 1:length(qs[i])

            out += qs[i][k](r, κs_λ[i])
        end
    end

    ## slower possibly due to r = u_rad - d[i] being evaluated every time qs is called.
    #out = sum( sum(qs[i][k](u_rad - d[i], κs_λ[i]) for k = 1:length(qs[i])) for i = 1:length(qs) )

    return out
end

function evalitpproxysys(qs::Vector{Vector{Function}},
    u_rad::T, x::SpinSysParamsType1{T}, κs_α::Vector{Vector{T}})::Complex{T} where T

    d = x.d
    κs_λ = x.κs_λ
    #κs_β = x.κs_β

    @assert length(d) == length(qs)

    out = zero(Complex{T})

    for i = 1:length(qs)
        r = u_rad - d[i]

        for k = 1:length(qs[i])

            out += κs_α[i][k]*qs[i][k](r, κs_λ[i])
        end
    end

    ## slower possibly due to r = u_rad - d[i] being evaluated every time qs is called.
    #out = sum( sum(qs[i][k](u_rad - d[i], κs_λ[i]) for k = 1:length(qs[i])) for i = 1:length(qs) )

    return out
end

function evalitpproxysys(qs::Vector{Vector{Function}},
    u_rad::T, x::SpinSysParamsType2{T})::Complex{T} where T

    d = x.d
    κs_λ = x.κs_λ
    #κs_β = x.κs_β

    @assert length(d) == length(qs)

    out = zero(Complex{T})

    #u_rad = 2*π*u
    for i = 1:length(qs)

        for k = 1:length(qs[i])
            r = u_rad - d[i][k]

            #out += qs[i][k](r, κs_λ[i][k], κs_β[i])
            out += qs[i][k](r, κs_λ[i][k])
        end
    end

    return out
end

function evalκitpproxysys(κ_α::Vector{Vector{T}}, qs::Vector{Vector{Function}},
    u_rad::T, x::SpinSysParamsType1{T})::Complex{T} where T

    d = x.d
    κs_λ = x.κs_λ
    #κs_β = x.κs_β

    @assert length(d) == length(qs)

    out = zero(Complex{T})

    #u_rad = 2*π*u
    for i = 1:length(qs)
        r = u_rad - d[i]

        for k = 1:length(qs[i])
            #out += κ_α[i][k]*qs[i][k](r, κs_λ[i], κs_β[i])
            out += κ_α[i][k]*qs[i][k](r, κs_λ[i])
        end
    end

    return out
end

function evalκitpproxysys(κ_α::Vector{Vector{T}}, qs::Vector{Vector{Function}},
    u_rad::T, x::SpinSysParamsType2{T})::Complex{T} where T

    d = x.d
    κs_λ = x.κs_λ

    @assert length(d) == length(qs)

    out = zero(Complex{T})

    #u_rad = 2*π*u
    for i = 1:length(qs)

        for k = 1:length(qs[i])
            r = u_rad - d[i][k]

            #out += κ_α[i][k]*qs[i][k](r, κs_λ[i][k], κs_β[i])
            out += κ_α[i][k]*qs[i][k](r, κs_λ[i][k])
        end
    end

    return out
end


###################### singlets.
function evalsinglets(u_rad::T, d::Vector{T}, αs_singlets::Vector{T}, Ωs_singlets,
    βs_singlets, λ0::T, λ_multipliers::Vector{T}) where T <: Real

    #u_rad = 2*π*u

    out = zero(Complex{T})
    for i = 1:length(αs_singlets)
        τ = u_rad - d[i]

        λ = λ0*λ_multipliers[i]
        Ω = Ωs_singlets[i]
        out += αs_singlets[i]*exp(im*βs_singlets[i])/(λ+im*(τ-Ω))
    end
    return out
end

function evalsinglets(u_rad::T, d::Vector{T},
    αs_singlets::Vector{T}, Ωs_singlets,
    βs_singlets, λ0::T, λ_multipliers::Vector{T},
    κ_α_singlets::Vector{T}) where T <: Real

    #u_rad = 2*π*u

    out = zero(Complex{T})
    for i = 1:length(αs_singlets)
        τ = u_rad - d[i]

        λ = λ0*λ_multipliers[i]
        Ω = Ωs_singlets[i]
        out += κ_α_singlets[i]*αs_singlets[i]*exp(im*βs_singlets[i])/(λ+im*(τ-Ω))
    end
    return out
end

###################### front end.
# with proxy.
function evalitpproxymixture(u_rad, As::Vector{SHType{T}},
    Bs::Vector{FIDModelType{T,SST}};
    w::Vector{T} = ones(T, length(As)))::Complex{T} where {T <: Real,SST}

    #u_rad = 2*π*u

    out = zero(Complex{T})

    for n = 1:length(As)
        out += w[n]*evalitpproxycompound(u_rad, As[n], Bs[n])
    end

    return out
end

# with proxy.
function evalitpproxycompound(u_rad, A::SHType{T}, B::FIDModelType{T,SST})::Complex{T} where {T <: Real, SST}

    #u_rad = 2*π*u

    out_sys = evalitpproxysys(B.qs, u_rad, B.ss_params)

    out_singlets = evalsinglets(u_rad, B.d_singlets, A.αs_singlets, A.Ωs_singlets,
    B.β_singlets, B.λ0, B.κs_λ_singlets)

    return out_sys + out_singlets
end


# with κ-proxy. refactor/remove this later.
function evalitpproxymixture(u_rad, As::Vector{SHType{T}},
    Es::Vector{καFIDModelType{T,SST}};
    w::Vector{T} = ones(T, length(As)))::Complex{T} where {T <: Real,SST}

    #u_rad = 2*π*u

    out = zero(Complex{T})

    for n = 1:length(Es)
        out += w[n]*evalitpproxycompound(u_rad, As[n], Es[n])
    end

    return out
end

# with κ-proxy.
function evalitpproxycompound(u_rad, A::SHType{T}, E::καFIDModelType{T,SST})::Complex{T} where {T <: Real, SST}

    #u_rad = 2*π*u

    out_sys = evalitpproxysys(E.core.qs, u_rad, E.core.ss_params, E.κs_α)

    out_singlets = evalsinglets(u_rad, E.core.d_singlets, A.αs_singlets, A.Ωs_singlets,
    E.core.β_singlets, E.core.λ0, E.core.κs_λ_singlets, E.κs_α_singlets)

    return out_sys + out_singlets
end

function findsimplecoherences(c2::Vector{T};
    atol = 1e-6)::Vector{Int} where T
    
    c2a = collect( abs.(c2[i]) for i = 1:length(c2))
    sum_c2a = collect(sum(c2a[i]) for i = 1:length(c2a))
    flags = isapprox.(sum_c2a, 1.0, atol = atol)

    return collect(1:length(c2))[flags]
end