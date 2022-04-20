
"""
getΩS(As::Vector{CompoundFIDType{T,SST}}) where {T,SST}

returns ΩS::Vector{Vector{Vector{T}}}
"""
function getΩS(As::Vector{CompoundFIDType{T,SST}}) where {T,SST}

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

"""
getPs( ΩS::Vector{Vector{Vector{T}}}, hz2ppmfunc) where T <: Real

returns Ps::Vector{Vector{Vector{T}}}
"""
function getPs( ΩS::Vector{Vector{Vector{T}}}, hz2ppmfunc) where T <: Real

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

#### for cost function.

# getPsnospininfo
function getPsnospininfo(As::Vector{NMRSpectraSimulator.CompoundFIDType{T,SST}}, hz2ppmfunc) where {T,SST}

    ΩS_ppm = Vector{Vector{T}}(undef, length(As))

    for (n,A) in enumerate(As)

        ΩS_ppm[n] = hz2ppmfunc.( combinevectors(A.Ωs) ./ (2*π) )

        tmp = hz2ppmfunc.( A.Ωs_singlets ./ (2*π) )
        push!(ΩS_ppm[n], tmp...)
    end

    return ΩS_ppm
end

function combinevectors(x::Vector{Vector{T}})::Vector{T} where T

    if isempty(x)
        return Vector{T}(undef, 0)
    end

    N = sum(length(x[i]) for i = 1:length(x))

    y = Vector{T}(undef,N)

    st_ind = 0
    fin_ind = 0
    for i = 1:length(x)
        st_ind = fin_ind + 1
        fin_ind = st_ind + length(x[i]) - 1

        y[st_ind:fin_ind] = x[i]
    end

    return y
end

function isnumericallyclose(x::T, y::T, tol = eps(T)*2) where T
    if abs(x-y) < tol
        return true
    end

    return false
end

"""
    convertcompactdomain(x::T, a::T, b::T, c::T, d::T)::T

converts compact domain x ∈ [a,b] to compact domain out ∈ [c,d].
"""
function convertcompactdomain(x::T, a::T, b::T, c::T, d::T)::T where T <: Real

    return (x-a)*(d-c)/(b-a)+c
end

function convertcompactdomain(x::Vector{T}, a::T, b::T, c::T, d::T)::Vector{T} where T <: Real

    return collect( convertcompactdomain(x[i], a, b, c, d) for i = 1:length(x) )
end


"""
prunelowsignalentries(x, threshold_factor::T, reduction_factor::Int)::Vector{Int} where T

keep the high signal portion and the reduced low signal portion.
keep_threshold = maximum(abs_x) * threshold_factor.
"""
function prunelowsignalentries(x, threshold_factor::T, reduction_factor::Int) where T
    
    abs_x = abs.(x)

    # find indices that have a low signal.
    keep_threshold = maximum(abs_x) * threshold_factor
    inds_for_reducing = findall(xx->xx<keep_threshold, abs_x)
    reduced_inds = inds_for_reducing[1:reduction_factor:length(inds_for_reducing)]

    # find indices that have a high signal.
    keep_inds = findall(xx->xx>keep_threshold, abs_x)

    # keep the high signal portion and the reduced low signal portion.
    inds = sort([keep_inds; reduced_inds])

    return inds, keep_inds, inds_for_reducing
end