################### spin system.

function evalitpproxysys(qs::Vector{Vector{Function}},
    u_rad::T, x::SpinSysFIDType1{T})::Complex{T} where T

    d = x.d
    κs_λ = x.κs_λ
    #κs_β = x.κs_β

    @assert length(d) == length(qs)

    out = zero(Complex{T})
    # out_real::T = zero(T) # slower.
    # out_imag::T = zero(T)

    #u_rad = 2*π*u
    for i = 1:length(qs)
        r = u_rad - d[i]

        for k = 1:length(qs[i])

            out += qs[i][k](r, κs_λ[i])
            # a, b = qs[i][k](r, κs_λ[i]) # slower.
            # out_real += a
            # out_imag += b
        end
    end
    #out = out_real + im*out_imag

    return out
end

function evalitpproxysys(qs::Vector{Vector{Function}},
    u_rad::T, x::SpinSysFIDType2{T})::Complex{T} where T

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
    u_rad::T, x::SpinSysFIDType1{T})::Complex{T} where T

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
    u_rad::T, x::SpinSysFIDType2{T})::Complex{T} where T

    d = x.d
    κs_λ = x.κs_λ
    #κs_β = x.κs_β

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

function evalκsinglets(u_rad::T, d::Vector{T},
    αs_singlets::Vector{T}, Ωs_singlets,
    βs_singlets, λ0::T, λ_multipliers::Vector{T};
    κ_α_singlets::Vector{T} = ones(T, length(αs_singlets))) where T <: Real

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
function evalitpproxymixture(u_rad, As::Vector{CompoundFIDType{T,SST}};
    w::Vector{T} = ones(T, length(As)))::Complex{T} where {T <: Real,SST}

    #u_rad = 2*π*u

    out = zero(Complex{T})

    for n = 1:length(As)
        out += w[n]*evalitpproxycompound(u_rad, As[n])
    end

    return out
end

# with proxy.
function evalitpproxycompound(u_rad, A::CompoundFIDType{T,SST})::Complex{T} where {T <: Real, SST}

    #u_rad = 2*π*u

    out_sys = evalitpproxysys(A.qs, u_rad, A.ss_params)

    out_singlets = evalsinglets(u_rad, A.d_singlets, A.αs_singlets, A.Ωs_singlets,
    A.β_singlets, A.λ0, A.κs_λ_singlets)

    return out_sys + out_singlets
end


# with κ-proxy. refactor/remove this later.
function evalitpproxymixture(u_rad, Es::Vector{κCompoundFIDType{T,SST}};
    w::Vector{T} = ones(T, length(Es)))::Complex{T} where {T <: Real, SST}

    #u_rad = 2*π*u

    out = zero(Complex{T})

    for n = 1:length(Es)
        out += w[n]*evalitpproxycompound(u_rad, Es[n])
    end

    return out
end

# with κ-proxy.
function evalitpproxycompound(u_rad, A::κCompoundFIDType{T,SST})::Complex{T} where {T <: Real, SST}

    #u_rad = 2*π*u

    out_sys = evalκitpproxysys(A.κs_α, A.core.qs, u_rad, A.core.ss_params)

    out_singlets = evalκsinglets(u_rad, A.core.d_singlets,
    A.core.αs_singlets, A.core.Ωs_singlets,
    A.core.β_singlets, A.core.λ0, A.core.κs_λ_singlets;
    κ_α_singlets = A.κs_α_singlets)

    return out_sys + out_singlets
end
