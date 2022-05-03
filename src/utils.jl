
"""
getΩS(As::Vector{CompoundFIDType{T,SST}}) where {T,SST}

returns ΩS::Vector{Vector{Vector{T}}}
"""
function getΩS(As::Vector{SHType{T}}) where T

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
function getPsnospininfo(As::Vector{NMRSpectraSimulator.SHType{T}}, hz2ppmfunc) where T

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


### for plotting.
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


### for reading config files.
function extractinfofromconfig( config_path::String,
    molecule_names::Vector{String}) where T <: Real

    N_compounds = length(molecule_names)

    # read.
    file_strings = readlines(config_path)

    cs_delta_group = Vector{Vector{Vector{Float64}}}(undef, N_compounds)
    #λ_group_labels = Vector{Vector{Vector{Int}}}(undef, N_compounds)

    for n = 1:N_compounds

        cs_delta_group[n] = extractmoleculeinfofromconfig(file_strings,
                                    molecule_names[n], config_path)
    end

    return cs_delta_group#, λ_group_labels
end

function extractmoleculeinfofromconfig( file_strings,
                                        name::String,
                                        file_label::String) where T <: Real

    # find the block of text for the input name.
    text_buffer = findblock(file_strings, name, "end compound", file_label)
    if isempty(text_buffer)
        println("Error processing config file.")
        return Vector{Vector{Float64}}(undef, 0), Vector{Vector{Int}}(undef, 0)
    end

    # get the number of groups.
    N_groups = length(filter(xx->occursin("Group",xx), text_buffer))

    # get default set up.
    cs_delta_group = Vector{Vector{Float64}}(undef, N_groups)
    #λ_group_labels = Vector{Vector{Int}}(undef, N_groups)

    for i = 1:N_groups

        # update position.
        cs_buffer = findblock(text_buffer, "Group $(i)", "end group", "Group $(i)")
        if isempty(cs_buffer)
            return Vector{Vector{Float64}}(undef, 0), Vector{Vector{Int}}(undef, 0)
        end
        N_cs = length(cs_buffer)

        cs_delta_group[i] = Vector{Float64}(undef, N_cs)
        #λ_group_labels[i] = Vector{Int}(undef, N_cs)

        for k = 1:length(cs_buffer)
            tokens = split(cs_buffer[k], (':',','))
            cs_delta_group[i][k] = tryparse(Float64, tokens[2])
            #λ_group_labels[i][k] = tryparse(Int, tokens[3])
        end
    end

    return cs_delta_group#, λ_group_labels
end


"""
uniqueinds(a_in::Vector{T}; atol = 1e-6) where T

Returns the unique values of `a_in`, with absolute tolerance `atol`, and the indices for each unique value.
"""
function uniqueinds(a_in::Vector{T}; atol = 1e-6) where T

    if length(a_in) < 2
        inds = Vector{Vector{Int}}(undef, length(a_in))
        
        if length(a_in) == 1
            inds[1] = Vector{Int}(undef, 1)
            inds[1][1] = 1
        end

        return copy(a_in), inds
    end

    a = copy(a_in)
    inds_a = collect(1:length(a))

    b = Vector{T}(undef, 0)
    inds_b = Vector{Vector{Int}}(undef, 0)

    while !isempty(a)

        target = pop!(a)
        inds = findall(xx->isapprox(target, xx; atol = atol), a)
        
        target_ind = pop!(inds_a)
        
        # save to output.
        push!(b, target)
        push!(inds_b, [target_ind; inds_a[inds]])

        # delete saved.
        deleteat!(a, inds)
        deleteat!(inds_a, inds)
    end

    return b, inds_b
end

"""
removetargetintegers(C::Vector{Vector{Int}}, search_list::Vector{Int})

remove an entry from C if the any of its integer array values appear in any entry of `search_list`.
"""
function removetargetintegers(C::Vector{Vector{Int}}, search_list::Vector{Int})

    keep_flags = trues(length(C))
    for i = 1:length(C)

        if any(C[i][1] .== search_list)
            keep_flags[i] = false
        end
    end


    return C[keep_flags]
end

"""
keeptargetintegers(C::Vector{Vector{Int}}, search_list::Vector{Int})

keep an entry from C if the any of its integer array values appear in any entry of `search_list`.
"""
function keeptargetintegers(C::Vector{Vector{Int}}, search_list::Vector{Int})

    keep_flags = falses(length(C))
    for i = 1:length(C)

        if any(C[i][1] .== search_list)
            keep_flags[i] = true
        end
    end


    return C[keep_flags]
end