# chemically equivalent nuclei are ones that have the same chemical shift.


### routines related to checking whether the J-coupling value between
#   chemically equivalent nuclei are the same.

"""
getpairs(inds::Vector{T})

get all exhaustive pairwise combos without symmetry of the 1D array `inds`.
"""
function getpairs(inds::Vector{T}) where T
    
    out = Vector{Tuple{T,T}}(undef, 0)
    for i = 1:length(inds)
        for j = i+1:length(inds)
            push!(out, (inds[i], inds[j]))
        end
    end

    return out
end

"""
getJfromdict(i::Int, j::Int, dict::Dict{Tuple{Int,Int},T} ) where T

query `dict` for the entry `(i,j)` and `(j,i)`. Returns the entry for `(i,j)`. Returns the entry for `(j,i)` if `(i,j)` is not found. Return zero if both entries are not found.
"""
function getJfromdict(i::Int, j::Int, dict::Dict{Tuple{Int,Int},T} ) where T
    #
    if haskey(dict, (i,j))
        return dict[(i,j)]
    end

    if haskey(dict, (j,i))
        return dict[(j,i)]
    end

    return zero(T)
end


"""
isallsame(a::Vector{T}; atol::T = 1e-6) where T

returns true if the entries in `a` are all within an abolute tolerance of `atol`.
"""
function isallsame(a::Vector{T}; atol::T = 1e-6) where T
    if length(findall(xx->isapprox(a[1], xx; atol = atol), a)) == length(a)
        return true
    end

    return false
end

#### misc.

"""
mapQtoHID(Q_g::Vector{Vector{Int}}, H_IDs::Vector{Int})

Convert entries in `Q_g` from taking on integer values between 1 to length(`H_IDs`) to taking on values in `H_IDs`.
"""
function mapQtoHID(Q_g::Vector{Vector{Int}}, H_IDs::Vector{Int})

    dict_g_to_H_IDs = Dict( collect(1:length(H_IDs)) .=> H_IDs)

    Q = Vector{Vector{Int}}(undef, length(Q_g))
    for i = 1:length(Q_g)

        Q[i] = Vector{Int}(undef, length(Q_g[i]))
        for k = 1:length(Q_g[i])
            Q[i][k] = dict_g_to_H_IDs[Q_g[i][k]]
        end
    end

    return Q
end

### routines related to checking whether the J-coupling value of A-C and B-C
# are the same, for all C's connected to either A or B. A and B are
# chemically equivalent nuclei.
"""
matchJlabels((label_pair_list::Vector{Tuple{Int,Int}}, search_list::Vector{Int})::Vector{Int}

Returns the indices of entries of `label_pair_list` that has any of the labels in the `search_list`.
"""
function matchanyJlabels(label_pair_list::Vector{Tuple{Int,Int}}, search_list::Vector{Int})::Vector{Int}

    inds = Vector{Int}(undef, 0)
    for k = 1:length(label_pair_list)
        i, j = label_pair_list[k]

        if any(i .== search_list) || any(j .== search_list)
            push!(inds, k)
        end
    end

    return inds
end

"""
getJIDstest(J_IDs::Vector{Tuple{Int,Int}},
    common_cs_IDs::Vector{Int},
    dict_H_IDs_to_css; atol = 1e-6)

Get the list of A-C, B-C, A-D, B-D, etc pairs to test.
"""
function getJIDstest(J_IDs::Vector{Tuple{Int,Int}},
    common_cs_IDs::Vector{Int},
    dict_H_IDs_to_css; atol = 1e-6)

    inds_any = matchanyJlabels(J_IDs, common_cs_IDs)
    inds_both = NMRSpectraSimulator.matchJlabels(J_IDs, common_cs_IDs)
    inds = setdiff(inds_any, inds_both)

    #inds_both = 
    J_IDs_tmp = J_IDs[inds]

    # remove the ID pairs that only common IDs that have the same chemical shift.
    J_IDs_test = Vector{Tuple{Int,Int}}(undef, 0)
    for l = 1:length(J_IDs_tmp)
        i, j = J_IDs_tmp[l]

        if !isapprox(dict_H_IDs_to_css[i], dict_H_IDs_to_css[j], atol = atol)
            push!(J_IDs_test, (i,j))
        end
    end

    return J_IDs_test
end

"""
checkmageq(test_inds::Vector{Vector{Int}},
    cs_IDs::Vector{Int},
    dict_J_ID_to_val;
    atol = 1e-6)::Bool

`test_inds` is a list of indices in cs_IDs. This function checks for magnetic equivalence for the nuclei in `test_inds`.
"""
function checkmageq(test_inds::Vector{Int},
    cs_IDs::Vector{Int},
    dict_J_ID_to_val;
    atol = 1e-6)::Bool

    if length(test_inds) == 1
        # need to have at least two nuclei to be magnetically equivalent.
        return false
    end

    
    common_cs_IDs = cs_IDs[test_inds]
    p_IDs = getpairs(common_cs_IDs)

    ## test common J within pairs.
    J_p = collect( getJfromdict(p_IDs[l][1], p_IDs[l][2], dict_J_ID_to_val) for l = 1:length(p_IDs) )
    pass_common_J_flag = isallsame(J_p; atol = atol)

    ## test common J with other IDs.
 
    # find all connections to first pair.
    t_IDs = getJIDstest(J_IDs, common_cs_IDs, dict_H_IDs_to_css; atol = atol)

    # check J-values.
    J_t = collect( getJfromdict(t_IDs[l][1], t_IDs[l][2], dict_J_ID_to_val) for l = 1:length(t_IDs) )
    pass_test_J_flag = isallsame(J_t; atol = atol)

    return pass_common_J_flag & pass_test_J_flag
end

"""
getmageqIDs(C::Vector{Vector{Int}},
    dict_ind_to_H_ID,
    dict_H_inds_to_css;
    cs_round_digits::Int = 5,
    atol::Float64 = 1e-6)

Given a list of node indices `C` of cliques of an undirected graph, with the indices of a clique being C[i],
Returns the list of node indices for each clique that are magnetically equivalent.
"""
function getmageqIDs(C::Vector{Vector{Int}},
    dict_ind_to_H_ID,
    dict_H_inds_to_css;
    cs_round_digits::Int = 5,
    atol::Float64 = 1e-6)

    # traverse each maximally connected cliques.
    C_IDs_mag_eq = Vector{Vector{Vector{Int}}}(undef, length(C))
    for i = 1:length(C)

        ### partition the nucleus by unique cs.
        cs = collect( dict_H_inds_to_css[C[i][k]] for k = 1:length(C[i]) )
        cs_IDs = collect( dict_ind_to_H_ID[C[i][k]] for k = 1:length(C[i]) )

        unique_cs = unique(round.(cs, digits = cs_round_digits))
        unique_cs_inds = collect( inds = findall(xx->isapprox(unique_cs[k], xx; atol = atol), cs) for k = 1:length(unique_cs) )

        # check magnetic equivalence for the chemically equivalent entries in `unique_cs_inds`.
        pass_flags = collect( checkmageq(unique_cs_inds[k], cs_IDs, dict_J_ID_to_val; atol = atol) for k = 1:length(unique_cs_inds) )

        # store result.
        common_cs_IDs = collect( cs_IDs[unique_cs_inds[k]] for k = 1:length(unique_cs_inds) )
        C_IDs_mag_eq[i] = common_cs_IDs[pass_flags]

        if length(common_cs_IDs[pass_flags]) > 2
            println("Warning for clique $(i): more than two groups of magnetically equivalent nuclei!")
        end
    end
    
    return C_IDs_mag_eq
end