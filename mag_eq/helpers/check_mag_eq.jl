# chemically equivalent nuclei are ones that have the same chemical shift.


### routines related to checking whether the J-coupling value between
#   chemically equivalent nuclei are the same.







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










