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
output index format:
x[i][j][k][l]
i-th spin system, k-th magnetically equivalent nuclei, l-th local spin index.
"""
function getmageqcompound(g,
    H_inds_sys,
    dict_ind_to_H_ID,
    dict_H_inds_to_css, 
    dict_H_IDs_to_css;
    atol = 1e-6)

    C_g = Graphs.maximal_cliques(g)

    N_spin_systems = length(H_inds_sys)
    
    mag_eq_sys_inds_local = Vector{Vector{Vector{Int}}}(undef, N_spin_systems)
    mag_eq_sys_inds_global = Vector{Vector{Vector{Int}}}(undef, N_spin_systems)
    mag_eq_sys_IDs = Vector{Vector{Vector{Int}}}(undef, N_spin_systems)

    for i = 1:N_spin_systems
        C = NMRSpectraSimulator.keeptargetintegers(C_g, H_inds_sys[i])

        mag_eq_IDs0, mag_eq_inds0 = getmageqIDs(C,
            dict_ind_to_H_ID,
            dict_H_inds_to_css;
            atol = atol)

        mag_eq_inds = combinetransitiveeqgroups(mag_eq_inds0)
        mag_eq_IDs = combinetransitiveeqgroups(mag_eq_IDs0)

        # mag_eq_inds = mag_eq_inds0
        # mag_eq_IDs = mag_eq_IDs0
        
        mag_eq_sys_inds_local[i] = getmageqlocalinds(H_inds_sys[i], mag_eq_inds)
        mag_eq_sys_IDs[i] = mag_eq_IDs
        mag_eq_sys_inds_global[i] = mag_eq_inds
    end

    return mag_eq_sys_inds_local, mag_eq_sys_IDs, mag_eq_sys_inds_global
end

"""
combinetransitiveeqgroups(eq_inds::Vector{Vector{Int}})

`eq_inds` contains theequivalent nuclei indices.
This function combine the groups of equivalent nuclei via transitivity.

Assumes the equivalence relation cannot be true for a group containg only one nucleus.
"""
function combinetransitiveeqgroups(eq_inds::Vector{Vector{Int}})

    if isempty(eq_inds)
        return Vector{Vector{Int}}(undef, 0)
    end

    eq_inds_flat_unique = unique(NMRSpectraSimulator.combinevectors(eq_inds))
    V = 1:length(eq_inds_flat_unique)
    dict_eq_to_V = Dict(eq_inds_flat_unique .=> V)
    dict_V_to_eq = Dict(V .=> eq_inds_flat_unique)
    
    eq_inds_V = collect( collect( dict_eq_to_V[eq_inds[i][k]] for k = 1:length(eq_inds[i])) for i = 1:length(eq_inds) )

    # get the edges of fully connected 
    E = Vector{Tuple{Int,Int}}(undef, 0)
    for k = 1:length(eq_inds)
        E_k = getconnectpath(eq_inds_V[k])

        push!(E, E_k...)
    end

    #
    h = Graphs.SimpleDiGraph(V[end])
    for k = 1:length(E)
        Graphs.add_edge!(h, E[k][1], E[k][2])
    end

    H = Graphs.connected_components(h)

    # convert back to the nuclei indexing of `eq_inds`.
    H_out = collect( collect( dict_V_to_eq[H[i][k]] for k = 1:length(H[i])) for i = 1:length(H) )

    return H_out
end

"""
Returns the edges of a connected path as specified by `vertices`.
"""
function getconnectpath(vertices::Vector{Int})
    N = length(vertices)

    out = Vector{Tuple{Int,Int}}(undef, N-1)
    
    for i = 1:N-1
        out[i] = (vertices[i], vertices[i+1])
    end

    return out
end

# """
# Returns the edges of the maximal clique specified by `vertices`.
# """
# function edgesofmaximalclique(vertices::Vector{Int})
#     N = length(vertices)

#     out = Vector{Tuple{Int,Int}}(undef, div(N*(N-1),2))
    
#     k = 0
#     for i = 1:N
#         for j = i+1:N

#             k += 1
#             out[k] = (vertices[i], vertices[j])
#         end
#     end

#     return out
# end


"""
getmageqIDs(C::Vector{Vector{Int}},
    dict_ind_to_H_ID,
    dict_H_inds_to_css;
    atol::Float64 = 1e-6)

Given a list of node indices `C` of cliques of an undirected graph, with the indices of a clique being C[i],
Returns the list of node indices for each clique that are magnetically equivalent.
"""
function getmageqIDs(C::Vector{Vector{Int}},
    dict_ind_to_H_ID,
    dict_H_inds_to_css;
    atol::Float64 = 1e-6)

    # traverse each maximally connected cliques.
    C_IDs_mag_eq = Vector{Vector{Vector{Int}}}(undef, length(C))
    C_inds_mag_eq = Vector{Vector{Vector{Int}}}(undef, length(C))
    for i = 1:length(C)

        ### partition the nucleus by unique cs.
        cs = collect( dict_H_inds_to_css[C[i][k]] for k = 1:length(C[i]) )
        cs_IDs = collect( dict_ind_to_H_ID[C[i][k]] for k = 1:length(C[i]) )
        cs_inds = C[i]

        # unique_cs = unique(round.(cs, digits = cs_round_digits))
        # unique_cs_inds = collect( inds = findall(xx->isapprox(unique_cs[k], xx; atol = atol), cs) for k = 1:length(unique_cs) )

        unique_cs, unique_cs_inds = NMRSpectraSimulator.uniqueinds(cs; atol = atol)

        # check magnetic equivalence for the chemically equivalent entries in `unique_cs_inds`.
        pass_flags = collect( checkmageq(unique_cs_inds[k], cs_IDs, dict_J_ID_to_val; atol = atol) for k = 1:length(unique_cs_inds) )

        # store result.
        common_cs_IDs = collect( cs_IDs[unique_cs_inds[k]] for k = 1:length(unique_cs_inds) )
        C_IDs_mag_eq[i] = common_cs_IDs[pass_flags]

        common_cs_inds = collect( cs_inds[unique_cs_inds[k]] for k = 1:length(unique_cs_inds) )
        C_inds_mag_eq[i] = common_cs_inds[pass_flags]

        if length(common_cs_IDs[pass_flags]) > 2
            println("Warning for clique $(i): more than two groups of magnetically equivalent nuclei!")
        end
    end

    # remove empty entries.
    keep_flags = trues(length(C))
    for i = 1:length(C)
        if isempty(C_IDs_mag_eq[i])
            keep_flags[i] = false
        end
    end
    C_IDs_mag_eq = C_IDs_mag_eq[keep_flags]
    C_inds_mag_eq = C_inds_mag_eq[keep_flags]

    mag_eq_IDs::Vector{Vector{Int}} = NMRSpectraSimulator.combinevectors(C_IDs_mag_eq)
    mag_eq_inds::Vector{Vector{Int}} = NMRSpectraSimulator.combinevectors(C_inds_mag_eq)
    
    return mag_eq_IDs, mag_eq_inds
end



"""
convertlabelsglobaltolocal(H_inds_sys::Vector{Vector{Int}},
    mag_eq_sys_inds::Vector{Vector{Int}})

The return type is Vector{Vector{Int}}.

Local: Each spin system will have its own spin nucleui numbering that starts at 1.
Global: The spins systems will together have one single spin nucleui numbering that starts at 1.
"""
function getmageqlocalinds(H_inds::Vector{Int},
    mag_eq_sys_inds::Vector{Vector{Int}})

    N_eqs = length(mag_eq_sys_inds)
    labels_sys_local = Vector{Vector{Int}}(undef, N_eqs)

    for i = 1:N_eqs

        x = H_inds
        conversion_dict = Dict(x .=> collect(1:length(x)))

        z = mag_eq_sys_inds[i]
        labels_sys_local[i] = collect( conversion_dict[z[k]] for k = 1:length(z) )  
    end

    return labels_sys_local
end


