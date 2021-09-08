
function getPs( ΩS::Vector{Vector{Vector{T}}},
                hz2ppmfunc) where T <: Real

    N_compounds = length(ΩS)

    Ps = Vector{Vector{Vector{T}}}(undef, N_compounds)
    for n = 1:N_compounds

        Ps[n] = Vector{Vector{T}}(undef, length(ΩS[n]))
        for i = 1:length(ΩS[n])

            Ps[n][i] = Vector{T}(undef, length(ΩS[n][i]))
            for l = 1:length(ΩS[n][i])

                Ps[n][i][l] = hz2ppmfunc( ΩS[n][i][l]/(2*π) )
            end
        end
    end

    return Ps
end

# get unique values of css_sys in graph.jl.
# It is denoted here as cs.
function cstopcs(   cs::Vector{Vector{T}},
                    cs_LUT::Vector{Vector{Vector{Int}}}) where T <: Real
    #
    N_subsystems = length(cs)
    p_cs = Vector{Vector{T}}(undef, N_subsystems)

    for m = 1:N_subsystems

        p_cs[m] = Vector{T}(undef, length(cs_LUT[m]))

        for i = 1:length(cs_LUT[m])
            inds = cs_LUT[m][i]

            # take average.
            p_cs[m][i] = sum( cs[m][k] for k in inds)/length(inds)
        end

    end

    return p_cs
end

function getcslengthfromLUT(cs_LUT_sys::Vector{Vector{Vector{Int}}})

    N_subsystems = length(cs_LUT_sys)

    cs_len_sys = Vector{Int}(undef, N_subsystems)

    for m = 1:N_subsystems
        cs_LUT = cs_LUT_sys[m]

        cs_len_sys[m] = maximum(maximum(cs_LUT[i]) for i = 1:length(cs_LUT))
    end

    return cs_len_sys
end

function pcstocs(   p_cs::Vector{Vector{T}},
                    cs_LUT::Vector{Vector{Vector{Int}}},
                    cs_len_sys::Vector{Int}) where T <: Real

    N_subsystems = length(p_cs)

    cs = Vector{Vector{T}}(undef, N_subsystems)
    pcstocs!(cs, p_cs, cs_LUT, cs_len_sys)

    return cs
end

function pcstocs!(  cs::Vector{Vector{T}},
                    p_cs::Vector{Vector{T}},
                    cs_LUT::Vector{Vector{Vector{Int}}},
                    cs_len_sys::Vector{Int}) where T <: Real
    #
    N_subsystems = length(p_cs)

    for m = 1:N_subsystems

        cs[m] = Vector{T}(undef, cs_len_sys[m])

        # for i = 1:length(p_cs[m])
        #     inds = cs_LUT[m][i]
        #
        #     for k in inds
        #         cs[m][k] = p_cs[m][i]
        #     end
        # end
        pcstocs!(cs[m], p_cs[m], cs_LUT[m])
    end

    return nothing
end

function pcstocs!(  cs_m::Vector{T},
                    p_cs_m::Vector{T},
                    cs_LUT_m::Vector{Vector{Int}}) where T <: Real
    #
    @assert length(p_cs_m) == length(cs_LUT_m) <= length(cs_m)

    for i = 1:length(p_cs_m)
        inds = cs_LUT_m[i]

        for k in inds
            cs_m[k] = p_cs_m[i]
        end
    end

    return nothing
end


# initialize FID parameter. For βS and λS.
function initializeFIDparameter(ΩS::Vector{Vector{Vector{T}}},
                                fill_val::T)::Vector{Vector{Vector{T}}} where T <: Real
    #
    N = length(ΩS)

    λS = Vector{Vector{Vector{T}}}(undef, N)
    for n = 1:N

        N_groups = length(ΩS[n])
        λS[n] = Vector{Vector{T}}(undef, N_groups)

        for i = 1:N_groups

            L = length(ΩS[n][i])
            λS[n][i] = Vector{T}(undef, L)
            fill!(λS[n][i], fill_val)
        end
    end

    return λS
end

function initializeFIDparameter(ΩS::Vector{Vector{Vector{T}}},
    lb::T, ub::T)::Vector{Vector{Vector{T}}} where T <: Real
    #
    N = length(ΩS)

    βS = Vector{Vector{Vector{T}}}(undef, N)
    for n = 1:N

        N_groups = length(ΩS[n])
        βS[n] = Vector{Vector{T}}(undef, N_groups)

        for i = 1:N_groups

            L = length(ΩS[n][i])
            βS[n][i] = Vector{T}(undef, L)
            for l = 1:L
                βS[n][i][l] = convertcompactdomain(rand(), 0.0, 1.0, lb, ub)
            end
        end
    end

    return βS
end
