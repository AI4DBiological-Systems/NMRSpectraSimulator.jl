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

    #real_itp = Interpolations.interpolate(real.(A), Interpolations.BSpline(Interpolations.Cubic(Interpolations.Line(Interpolations.OnGrid()))))
    real_itp = Interpolations.interpolate(real.(A), Interpolations.BSpline(Interpolations.Quadratic(Interpolations.Line(Interpolations.OnGrid()))))
    real_sitp = Interpolations.scale(real_itp, A_r, A_ξ)
    real_setp = Interpolations.extrapolate(real_sitp, 0.0) # zero outside interp range.

    #imag_itp = Interpolations.interpolate(imag.(A), Interpolations.BSpline(Interpolations.Cubic(Interpolations.Line(Interpolations.OnGrid()))))
    imag_itp = Interpolations.interpolate(imag.(A), Interpolations.BSpline(Interpolations.Quadratic(Interpolations.Line(Interpolations.OnGrid()))))
    imag_sitp = Interpolations.scale(imag_itp, A_r, A_ξ)
    imag_setp = Interpolations.extrapolate(imag_sitp, 0.0) # zero outside interp range.

    #return real_sitp, imag_sitp
    return real_setp, imag_setp
end

function setupcompoundpartitionitp(d_max::T,
    x::SST,
    Δc_m_compound::Vector{Vector{Vector{T}}},
    part_inds_compound::Vector{Vector{Vector{Int}}},
    αs::Vector{Vector{T}}, Ωs::Vector{Vector{T}},
    λ0::T, u_min::T, u_max::T;
    κ_λ_lb = 0.5,
    κ_λ_ub = 2.5,
    Δr = 1.0,
    Δκ_λ = 0.05) where {T <: Real, SST}

    qs = Vector{Vector{Function}}(undef, length(αs))
    gs = Vector{Vector{Function}}(undef, length(αs)) # no phase.
    Δc_avg = Vector{Vector{Vector{T}}}(undef, length(αs))

    for i = 1:length(αs) # over elements in a spin group.

        N_partition_elements = length(part_inds_compound[i])
        qs[i] = Vector{Function}(undef, N_partition_elements)
        gs[i] = Vector{Function}(undef, N_partition_elements)
        Δc_avg[i] = Vector{Vector{T}}(undef, N_partition_elements)

        for k = 1:N_partition_elements
            #println("i,k", (i,k))

            inds = part_inds_compound[i][k]
            α = αs[i][inds]
            Ω = Ωs[i][inds]
            real_sitp, imag_sitp = setuppartitionitp(α, Ω,
            d_max, λ0, u_min, u_max; κ_λ_lb = κ_λ_lb, κ_λ_ub = κ_λ_ub,
            Δr = Δr, Δκ_λ = Δκ_λ)

            Δc_avg[i][k] = Statistics.mean( Δc_m_compound[i][inds] )
            #qs[i][k] = (rr, ξξ, bb)->(real_sitp(rr,ξξ)+im*imag_sitp(rr,ξξ))*exp(im*dot(bb, c)) # unpackaged.
            #qs[i][k] = (rr, ξξ, bb)->evalq(real_sitp, imag_sitp, rr, ξξ, bb, Δc_avg[i][k]) # packaged.
            qs[i][k] = (rr, ξξ)->evalq(real_sitp, imag_sitp, rr, ξξ, x.κs_β[i], Δc_avg[i][k])

            # x_buffer = Vector{T}(undef, 1)
            # s_x_buffer = Vector{T}(undef, 1)
            # c_x_buffer = Vector{T}(undef, 1)
            # a_buffer = Vector{T}(undef, 1)
            # b_buffer = Vector{T}(undef, 1)
            # qs[i][k] = (rr, ξξ)->evalq1!(x_buffer, s_x_buffer, c_x_buffer, a_buffer, b_buffer,
            #     real_sitp, imag_sitp, rr, ξξ, x.κs_β[i], Δc_avg[i][k])

            #qs[i][k] = (rr, ξξ, bb)->(real_sitp(rr,ξξ)+im*imag_sitp(rr,ξξ))*exp(im*dot(bb, Δc_m_compound[i][k]))
            #qs[i][k] = (rr, ξξ, bb)->(real_sitp(rr,ξξ)+im*imag_sitp(rr,ξξ))*exp(im*bb)

            gs[i][k] = (rr, ξξ)->(real_sitp(rr,ξξ)+im*imag_sitp(rr,ξξ))

            #gs[i][k] = (rr, ξξ)->real_sitp(rr,ξξ) # timing test for itp.
            #gs[i][k] = (rr, ξξ)->imag_sitp(rr,ξξ) # timing test for itp.

            # qrs[i][k] = real_sitp
            # qis[i][k] = imag_sitp
        end
    end

    return qs, gs, Δc_avg
end

function evalq(real_sitp, imag_sitp, r::T, ξ::T, b::Vector{T}, c)::Complex{T} where T <: Real

    return (real_sitp(r,ξ)+im*imag_sitp(r,ξ))*exp(im*dot(b, c))

    # # i (a sin(x) + b cos(x)) + a cos(x) - b sin(x)
    # x = dot(b,c)
    # s_x = sin(x)
    # c_x = cos(x)
    # a = real_sitp(r,ξ)
    # b = imag_sitp(r,ξ)
    # return im*(a*s_x + b*c_x) + a*c_x - b*s_x

    ## speed test.
    #return (real_sitp(r,ξ)+im*imag_sitp(r,ξ))
    #return real_sitp(r,ξ)
end

# # actually slower.
# function evalq1!(x::Vector{T},
#     s_x::Vector{T},
#     c_x::Vector{T},
#     a::Vector{T},
#     b::Vector{T},
#     real_sitp, imag_sitp, r::T, ξ::T, p_β::Vector{T}, c)::Tuple{T,T} where T <: Real
#
#     #return (real_sitp(r,ξ)+im*imag_sitp(r,ξ))*exp(im*dot(b, c))
#
#     # i (a sin(x) + b cos(x)) + a cos(x) - b sin(x)
#     x[1] = dot(p_β, c)
#     s_x[1] = sin(x[1])
#     c_x[1] = cos(x[1])
#     a[1] = real_sitp(r, ξ)
#     b[1] = imag_sitp(r, ξ)
#     #return im*(a[1]*s_x[1] + b[1]*c_x[1]) + a[1]*c_x[1] - b[1]*s_x[1]
#     return a[1]*c_x[1] - b[1]*s_x[1], a[1]*s_x[1] + b[1]*c_x[1]
#
#     ## speed test.
#     #return (real_sitp(r,ξ)+im*imag_sitp(r,ξ))
#     #return real_sitp(r,ξ)
# end
