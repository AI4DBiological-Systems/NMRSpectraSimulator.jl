# for 1D. The order of points in X matters.
function fuseptsgetID(   X_in::Vector{T};
                    radius = 0.1,
                    Minkowski_parameter = 3.5) where T <: Real
    X = X_in'
    N = length(X_in)

    # set up tree for fast search.
    #kdtree = NearestNeighbors.KDTree(X; leafsize = 10)
    balltree = NearestNeighbors.BallTree(X,
                    NearestNeighbors.Minkowski(Minkowski_parameter);
                    reorder = false)

    # initialize the group for each entry in X.
    group_IDs = zeros(Int, N)

    j = 0
    for i = 1:N

        if group_IDs[i] < 1
            # X[i] is an unprocessed point.

            inds = NearestNeighbors.inrange(balltree, [X[i]], radius, true)

            if !isempty(inds)
                # new group ID.
                j += 1
            end

            for l = 1:length(inds)


                group_IDs[inds[l]] = j

            end
        end
    end
    L = j # number of unique IDs in group_IDs.

    return group_IDs, L
end

function fuseΩα( Ωs::Vector{T},
                αs::Vector{T},
                group_IDs::Vector{Int},
                k::Int)::Tuple{T, Int, T} where T
    N = length(αs)
    @assert length(Ωs) == length(group_IDs) == N

    counter = 0
    numerator = zero(T)
    sum_α = zero(T)

    for i = 1:N

        if group_IDs[i] == k
            counter += 1

            # # equally weighted.
            # numerator += Ωs[i]

            # weighted by αs[i]
            numerator += Ωs[i]*αs[i]

            sum_α += αs[i]
        end
    end

    # # equally weighted.
    # fused_Ω = numerator/counter

    # weighted by αs
    fused_Ω = numerator/sum_α

    return fused_Ω, counter, sum_α
end

# Assumes that the first entry in each tuple is the smaller one.
function getΩfroms(s::Vector{T},
                    state_labels::Vector{Tuple{Int,Int}}) where T
    #
    L = length(state_labels)

    F = Vector{T}(undef, L)
    for l = 1:L
        q, r = state_labels[l]
        F[l] = s[r] - s[q]
    end

    return F
end

# approximate Ωs and αs with less entries, by fusing entries in pts
#   within a specified radius.
# pts is a version of Ωs where the distance is used by the fusing operation.
function approxΩsαs(    pts,
                        Ωs::Vector{T},
                        αs::Vector{T};
                        radius = 0.1,
                        Minkowski_parameter = 3.5) where T <: Real
    #
    @assert length(pts) == length(Ωs) == length(αs)

    group_IDs, L = fuseptsgetID(pts;
                            radius = radius,
                            Minkowski_parameter = Minkowski_parameter)

    # find average of each group
    F = Vector{T}(undef, L)
    fill!(F, Inf) # debug.

    a = Vector{T}(undef, L)
    fill!(a, Inf) # debug.

    status_flags = falses(L)

    for l = 1:L
        fused_Ω, N_l, sum_α = fuseΩα(Ωs, αs, group_IDs, l)

        # some group_IDs get overwritten,
        #   so need to check if this l is still valid.
        if N_l > 0
            F[l] = fused_Ω
            a[l] = sum_α
            status_flags[l] = true
        end
    end
    F = F[status_flags]
    a = a[status_flags]

    return F, a
end


# for each l ∈ [L],
#   - find the closest entry in F that approx. each Ω[l].
#   - return its corresponding entry in state_labels_F
function findapproxstates(  F::Vector{T},
                            state_labels_F::Vector{Tuple{Int,Int}},
                            Ω::Vector{T}) where T

    #
    L = length(Ω)
    N = length(F)

    state_labels = Vector{Tuple{Int,Int}}(undef, L)

    for l = 1:L
        dists = collect( norm(F[n]-Ω[l]) for n = 1:N)
        _, ind = findmin(dists)

        state_labels[l] = state_labels_F[ind]
    end

    return state_labels
end
