



function evalcompound(u, A::CompoundFIDType{T})::Complex{T} where T <: Real

    u_rad = 2*π*u

    out_sys = evalcLcompoundviapartitions(u, A.d, A.αs, A.Ωs, A.κs_λ, A.κs_β,
    A.λ0, A.Δc_m_compound, A.part_inds_compound)

    out_singlets = evalsinglets(u, A.d_singlets, A.αs_singlets, A.Ωs_singlets,
    A.β_singlets, A.λ0, A.κs_λ_singlets)

    return out_sys + out_singlets
end

function evalmixture(u, As::Vector{CompoundFIDType{T}};
    w::Vector{T} = ones(T, length(As)))::Complex{T} where T <: Real

    u_rad = 2*π*u

    out = zero(Complex{T})
    for n = 1:length(As)
        out += w[n]*evalcompound(u, As[n])
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

function evalcLcompoundviapartitions(u, d::Vector{T},
    αs::Vector{Vector{T}}, Ωs::Vector{Vector{T}},
    κs_λ::Vector{T}, κs_β::Vector{Vector{T}}, λ0::T,
    Δc_m_compound, part_inds_compound)::Complex{T} where T <: Real

    u_rad = 2*π*u

    out = zero(Complex{T})
    for i = 1:length(αs)
        r = u_rad-d[i]

        for k = 1:length(part_inds_compound[i])
            #r = 2*π*u-d[i][k]
            inds = part_inds_compound[i][k]

            out += evalcLpartitionelement(r, αs[i][inds],
            Ωs[i][inds], κs_λ[i]*λ0)*exp(im*dot(κs_β[i], Δc_m_compound[i][k]))
        end
    end

    return out
end


function setuppartitionitp(α::Vector{T}, Ω::Vector{T}, d_max::T, λ0::T,
    u_min::T, u_max::T;
    κ_λ_lb = 0.5,
    κ_λ_ub = 2.5,
    Δr = 1.0,
    Δκ_λ = 0.05) where T <: Real

    # A_x1 = 1:.1:10
    # A_x2 = 1:.5:20
    # f(x1, x2) = log(x1+x2)
    # A = [f(x1,x2) for x1 in A_x1, x2 in A_x2]
    # itp = Interpolations.interpolate(A, BSpline(Cubic(Line(OnGrid()))))
    # sitp = Interpolations.scale(itp, A_x1, A_x2)
    # sitp(5., 10.) # exactly log(5 + 10)
    # sitp(5.6, 7.1) # approximately log(5.6 + 7.1)

    # set up bounds.

    r_min = 2*π*u_min - d_max
    r_max = 2*π*u_max + d_max

    # set up samples.
    A_r = r_min:Δr:r_max # large range.
    #A_r = r_min:20.0:r_max # large range.
    A_ξ = κ_λ_lb:Δκ_λ:κ_λ_ub
    # println("length(A_r) = ", length(A_r))
    # println("length(A_ξ) = ", length(A_ξ))

    # complex.
    f = (rr,ξξ)->evalcLpartitionelement(rr, α, Ω, ξξ*λ0)
    A = [f(x1,x2) for x1 in A_r, x2 in A_ξ]

    real_itp = Interpolations.interpolate(real.(A), Interpolations.BSpline(Interpolations.Cubic(Interpolations.Line(Interpolations.OnGrid()))))
    real_sitp = Interpolations.scale(real_itp, A_r, A_ξ)
    real_setp = Interpolations.extrapolate(real_sitp, 0.0) # zero outside interp range.

    imag_itp = Interpolations.interpolate(imag.(A), Interpolations.BSpline(Interpolations.Cubic(Interpolations.Line(Interpolations.OnGrid()))))
    imag_sitp = Interpolations.scale(imag_itp, A_r, A_ξ)
    imag_setp = Interpolations.extrapolate(imag_sitp, 0.0) # zero outside interp range.

    #return real_sitp, imag_sitp
    return real_setp, imag_setp
end

function setupcompoundpartitionitp(d_max::T,
    Δc_m_compound::Vector{Vector{Vector{T}}},
    part_inds_compound::Vector{Vector{Vector{Int}}},
    αs::Vector{Vector{T}}, Ωs::Vector{Vector{T}},
    λ0::T, u_min::T, u_max::T;
    κ_λ_lb = 0.5,
    κ_λ_ub = 2.5,
    Δr = 1.0,
    Δκ_λ = 0.05) where T <: Real

    qs = Vector{Vector{Function}}(undef, length(αs))
    for i = 1:length(αs) # over elements in a spin group.

        N_partition_elements = length(part_inds_compound[i])
        qs[i] = Vector{Function}(undef, N_partition_elements)

        for k = 1:N_partition_elements
            #println("i,k", (i,k))
            α = αs[i][part_inds_compound[i][k]]
            Ω = Ωs[i][part_inds_compound[i][k]]
            real_sitp, imag_sitp = setuppartitionitp(α, Ω,
            d_max, λ0, u_min, u_max; κ_λ_lb = κ_λ_lb, κ_λ_ub = κ_λ_ub,
            Δr = Δr, Δκ_λ = 0.05)

            qs[i][k] = (rr, ξξ, bb)->(real_sitp(rr,ξξ)+im*imag_sitp(rr,ξξ))*exp(im*dot(bb, Δc_m_compound[i][k]))
            #qs[i][k] = (rr, ξξ, bb)->(real_sitp(rr,ξξ)+im*imag_sitp(rr,ξξ))*exp(im*bb)

        end
    end

    return qs
end

function evalsinglets(u::T, d::Vector{T}, αs_singlets::Vector{T}, Ωs_singlets,
    βs_singlets, λ0::T, λ_multipliers::Vector{T}) where T <: Real

    u_rad = 2*π*u

    out = zero(Complex{T})
    for i = 1:length(αs_singlets)
        τ = u_rad - d[i]

        λ = λ0*λ_multipliers[i]
        Ω = Ωs_singlets[i]
        out += αs_singlets[i]*exp(im*βs_singlets[i])/(λ+im*(τ-Ω))
    end
    return out
end

function evalκsinglets(u::T, d::Vector{T},
    αs_singlets::Vector{T}, Ωs_singlets,
    βs_singlets, λ0::T, λ_multipliers::Vector{T};
    κ_α_singlets::Vector{T} = ones(T, αs_singlets)) where T <: Real

    u_rad = 2*π*u

    out = zero(Complex{T})
    for i = 1:length(αs_singlets)
        τ = u_rad - d[i]

        λ = λ0*λ_multipliers[i]
        Ω = Ωs_singlets[i]
        out += κ_α_singlets[i]*αs_singlets[i]*exp(im*βs_singlets[i])/(λ+im*(τ-Ω))
    end
    return out
end

function evalitpproxysys(qs::Vector{Vector{Function}},
    u::T, d::Vector{T}, κs_λ::Vector{T}, κs_β::Vector{Vector{T}})::Complex{T} where T

    @assert length(d) == length(qs)

    out = zero(Complex{T})

    u_rad = 2*π*u
    for i = 1:length(qs)
        r = u_rad - d[i]

        for k = 1:length(qs[i])
            out += qs[i][k](r, κs_λ[i], κs_β[i])
        end
    end

    return out
end

function evalκitpproxysys(κ_α::Vector{Vector{T}}, qs::Vector{Vector{Function}},
    u::T, d::Vector{T}, κs_λ::Vector{T}, κs_β::Vector{Vector{T}})::Complex{T} where T

    @assert length(d) == length(qs)

    out = zero(Complex{T})

    u_rad = 2*π*u
    for i = 1:length(qs)
        r = u_rad - d[i]

        for k = 1:length(qs[i])
            out += κ_α[i][k]*qs[i][k](r, κs_λ[i], κs_β[i])
        end
    end

    return out
end

# with proxy.
function evalitpproxymixture(u, As::Vector{CompoundFIDType{T}};
    w::Vector{T} = ones(T, length(As)))::Complex{T} where T <: Real

    u_rad = 2*π*u

    out = zero(Complex{T})

    for n = 1:length(As)
        out += w[n]*evalitpproxycompound(u, As[n])
    end

    return out
end

# with proxy.
function evalitpproxycompound(u, A::CompoundFIDType{T})::Complex{T} where T <: Real

    u_rad = 2*π*u

    out_sys = evalitpproxysys(A.qs, u, A.d, A.κs_λ, A.κs_β)

    out_singlets = evalsinglets(u, A.d_singlets, A.αs_singlets, A.Ωs_singlets,
    A.β_singlets, A.λ0, A.κs_λ_singlets)

    return out_sys + out_singlets
end


#####

# with κ-proxy. refactor/remove this later.
function evalitpproxymixture(u, As::Vector{κCompoundFIDType{T}};
    w::Vector{T} = ones(T, length(As)))::Complex{T} where T <: Real

    u_rad = 2*π*u

    out = zero(Complex{T})

    for n = 1:length(As)
        out += w[n]*evalitpproxycompound(u, As[n])
    end

    return out
end

# with κ-proxy.
function evalitpproxycompound(u, A::κCompoundFIDType{T})::Complex{T} where T <: Real

    u_rad = 2*π*u

    out_sys = evalκitpproxysys(A.κ, A.core.qs, u, A.core.d, A.core.κs_λ, A.core.κs_β)

    out_singlets = evalκsinglets(u, A.core.d_singlets,
    A.core.αs_singlets, A.core.Ωs_singlets,
    A.core.β_singlets, A.core.λ0, A.core.κs_λ_singlets;
    κ_singlets = A.κ_singlets)

    return out_sys + out_singlets
end
