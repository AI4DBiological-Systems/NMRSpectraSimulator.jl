function getIja(Ia::Matrix{T}, j::Int, N::Int) where T

    Id = [one(T) zero(T); zero(T) one(T)]

    Ias = collect( Id for i = 1:N )
    Ias[j] = Ia
    return Kronecker.kron(Ias...)
end

function productIaIb(Ixs, Iys, Izs, a::Int, b::Int)
    return Ixs[a]*Ixs[b] + Iys[a]*Iys[b] + Izs[a]*Izs[b]
end

"""
For operator A.

Example:
N = 3
Ixs = collect( getIja(Ix, i, N) for i = 1:N )
Iys = collect( getIja(Iy, i, N) for i = 1:N )
Izs = collect( getIja(Iz, i, N) for i = 1:N )

# sanity check.
out1 = productIaIb(Ixs, Iys, Izs, 1, 3)
out2 = productIjIk(Ixs, Iys, Izs, [1;], [3;])
norm(out1-out2)
"""
function productIjIk(Ixs, Iys, Izs, inds_j::Vector{Int}, inds_k::Vector{Int})

    #

    out = sum(Ixs[j] for j in inds_j) * sum(Ixs[k] for k in inds_k) +
        sum(Iys[j] for j in inds_j) * sum(Iys[k] for k in inds_k) +
        sum(Izs[j] for j in inds_j) * sum(Izs[k] for k in inds_k)

    return out
end


"""
for operator B.

Example:
N = 3
Ixs = collect( getIja(Ix, i, N) for i = 1:N )
Iys = collect( getIja(Iy, i, N) for i = 1:N )
Izs = collect( getIja(Iz, i, N) for i = 1:N )

# # sanity check.
# out1 = productIaIb(Ixs, Iys, Izs, 1, 3)
# out2 = productIjIk(Ixs, Iys, Izs, [1;], [3;])

out1 = pairwiseproductIa(Ixs, Iys, Izs, [1; 2; 3])
out2 = productIaIb(Ixs, Iys, Izs, 1, 3) + productIaIb(Ixs, Iys, Izs, 2, 3) + productIaIb(Ixs, Iys, Izs, 1, 2)
norm(out1-out2)
"""
function pairwiseproductIa(Ixs, Iys, Izs, inds::Vector{Int})

    out = zeros(Float64, size(Ixs[1]))
    for i = 1:length(inds)

        ind1 = inds[i]
        for j = i+1:length(inds)
            ind2 = inds[j]

            out += productIaIb(Ixs, Iys, Izs, ind1, ind2)
        end
    end

    return out
end
