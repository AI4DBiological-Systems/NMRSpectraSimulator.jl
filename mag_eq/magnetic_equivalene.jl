using LinearAlgebra

include("../src/GISSMOReader.jl")
import .GISSMOReader

import NMRSpectraSimulator

import Kronecker
import Graphs

include("./helpers/graphs.jl")
include("./helpers/operators.jl")


# machine settings.
fs = 14005.602240896402
SW = 20.0041938620844
ν_0ppm = 10656.011933076665

hz2ppmfunc = uu->(uu - ν_0ppm)*SW/fs
ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)


# operators in Zeeman basis.
Id = NMRSpectraSimulator.getsingleId()
Ix = NMRSpectraSimulator.getsingleIx()
Iy_no_im = NMRSpectraSimulator.getsingleIynoim()
Iy = NMRSpectraSimulator.getsingleIy()
Iz = NMRSpectraSimulator.getsingleIz()




N = 4
Ixs = collect( getIja(Ix, i, N) for i = 1:N )
Iys = collect( getIja(Iy, i, N) for i = 1:N )
Izs = collect( getIja(Iz, i, N) for i = 1:N )

# # sanity check.
# out1 = productIaIb(Ixs, Iys, Izs, 1, 3)
# out2 = productIjIk(Ixs, Iys, Izs, [1;], [3;])

# # sanity check.
# out1 = pairwiseproductIa(Ixs, Iys, Izs, [1; 2; 3])
# out2 = productIaIb(Ixs, Iys, Izs, 1, 3) + productIaIb(Ixs, Iys, Izs, 2, 3) + productIaIb(Ixs, Iys, Izs, 1, 2)

#inds = [ [1;2;3]; [4;]; ]
#css = [3.4; 4.1]

# user input.
H_IDs = collect(1:4)
H_css = [2.3; 2.3; 2.3; 4.1; 4.1]

J_IDs = [ (1,2); (1,3); (1,4); (2,3); (2,4); (3,4)]
J_vals = [ 3.4; 3.4; 3.4; -5.8; -5.8; -5.8]

#load_path = "/home/roy/Documents/repo/GISSMOReader.jl/coupling_info/Maybridge_Ro3_Fragment_10_A08_simulation_1.json"
load_path = "/home/roy/Documents/repo/GISSMOReader.jl/coupling_info/bmse000860_simulation_1.json"
H_IDs, H_css, J_IDs, J_vals = GISSMOReader.loadcouplinginfojson(load_path)


# end user input.

ω_css = ppm2hzfunc.(H_css)

H_inds = collect(1:length(H_IDs))
J_inds = convertJIDstoJinds(J_IDs, H_IDs)

dict_H_ID_to_ind = Dict(H_IDs .=> collect(1:length(H_IDs)))
dict_ind_to_H_ID = Dict( collect(1:length(H_IDs)) .=> H_IDs)

dict_J_ID_to_val = Dict(J_IDs .=> J_vals)
dict_J_ind_to_val = Dict(J_inds .=> J_vals)

dict_J_ID_to_ind = Dict(J_IDs .=> J_inds)
dict_J_ind_to_ID = Dict(J_inds .=> J_IDs)


g = constructspinsystemsgraphforcompound(length(H_IDs), J_inds)

systems_g = Graphs.connected_components(g) # node inds for spin systems.
C = Graphs.maximal_cliques(g)

nodes = C[1]
# I am here. partition the nodes based on J_values.

# A1 = sum(w[k] .* sum(Izs[inds[k]]) for k = 1:length(w))
# A2 =
# B = Ixs[1]*Ixs[2]
# A36 = (A*B - B*A) # equation A.36 of Levitt's Spin Dynamics, "J-couplings and magnetic equivalence".
# println("N = $(N), norm(A36) = ", norm(A36))

@assert 1==2

N = 3
Ixs = collect( getIja(Ix, i, N) for i = 1:N )
Iys = collect( getIja(Iy, i, N) for i = 1:N )
Izs = collect( getIja(Iz, i, N) for i = 1:N )


A = (Ixs[1]+Ixs[2])*Ixs[3]
B = Ixs[1]*Ixs[2]
A36 = (A*B - B*A) # equation A.36 of Levitt's Spin Dynamics, "J-couplings and magnetic equivalence".
println("N = $(N), norm(A36) = ", norm(A36))

## book example of magnetic equivalence.
w12 = 11474.887839493498
w3 = 13205.050246451618
J13 = 5.8
J23 = 5.8
J12 = -7.9

# ## explore.
# w12 = 11474.887839493498
# w3 = 13205.050246451618
# #w3 = w12
# J13 = 0.0
# J23 = -7.9
# J12 = 0.0

A = w12*(Izs[1]+Izs[2]) + w3*Izs[3] + 2*π*J13*(Ixs[1]*Ixs[3] + Iys[1]*Iys[3] + Izs[1]*Izs[3]) +
    2*π*J23*(Ixs[2]*Ixs[3] + Iys[2]*Iys[3] + Izs[2]*Izs[3])
#B = real.(2*π*J12*(Ixs[1]*Ixs[2] + Iys[1]*Iys[2] + Izs[1]*Izs[2]))
B = 2*π*J12*(Ixs[1]*Ixs[2] + Iys[1]*Iys[2] + Izs[1]*Izs[2])

result = A*B - B*A # claim by Levitt just before equation A.36 of A and B commuting.
println("N = $(N), norm([A,B]) = ", norm(result))

Q = sum(Ixs)
result = B*Q - Q*B
println("N = $(N), norm([B,Q]) = ", norm(result))
println()

# observation operators (total angular moment along an axis) is always commutative to B.
C = Ixs[1]+Ixs[2]+Ixs[3]
Q_x = B*C - C*B # always 0.

C = Iys[1]+Iys[2]+Iys[3]
Q_y = B*C - C*B # always 0.

C = Izs[1]+Izs[2]+Izs[3]
Q_z = B*C - C*B # always 0.


Q = sum(Ixs)
#B = 2*π*J12*(pairwiseproductIa(Ixs, Iys, Izs, [1;2;]))
B2 = 2*π*J12*(productIaIb(Ixs, Iys, Izs, 1, 2)) #+ 2*π*(-5.2)*(productIaIb(Ixs, Iys, Izs, 2, 3))
result = B2*Q - Q*B2
println("N = $(N), norm([B2,Q]) = ", norm(result))

result = A*B2 - B2*A
println("N = $(N), norm([A,B2]) = ", norm(result))
println()

@assert 1==43

### book example, 4 spins.

N = 4
Ixs = collect( getIja(Ix, i, N) for i = 1:N )
Iys = collect( getIja(Iy, i, N) for i = 1:N )
Izs = collect( getIja(Iz, i, N) for i = 1:N )

A = sum(Ixs[1:3])*Ixs[4]
B = Ixs[1]*Ixs[2] + Ixs[2]*Ixs[3] + Ixs[1]*Ixs[3]
A36 = (A*B - B*A) # equation A.36 of Levitt's Spin Dynamics, "J-couplings and magnetic equivalence".
println("N = $(N), norm(A36) = ", norm(A36))

# book example of magnetic equivalence.
w123 = 11474.887839493498
w4 = 13205.050246451618
J12 = -7.9
J13 = J12
J23 = J12

J14 = 5.8
J24 = J14
J34 = J14



A = w123*sum(Izs[1:3]) + w4*Izs[4] + 2*π*J12*(addIa(Ixs, Iys, Izs, [1;2;3]), )
#B = real.(2*π*J12*(Ixs[1]*Ixs[2] + Iys[1]*Iys[2] + Izs[1]*Izs[2]))
B = 2*π*J12*(Ixs[1]*Ixs[2] + Iys[1]*Iys[2] + Izs[1]*Izs[2])

result = A*B - B*A # claim by Levitt just before equation A.36 of A and B commuting.
println("N = $(N), norm(result) = ", norm(result))

# observation operators (total angular moment along an axis) is always commutative to B.
C = Ixs[1]+Ixs[2]+Ixs[3]
Q_x = B*C - C*B # always 0.

C = Iys[1]+Iys[2]+Iys[3]
Q_y = B*C - C*B # always 0.

C = Izs[1]+Izs[2]+Izs[3]
Q_z = B*C - C*B # always 0.


### arbitrary number of spins, common J coupling values.


@assert 1==2444


# test loading.
load_path = "/home/roy/MEGAsync/inputs/NMR/debug/ethanol.json"

H_IDs, H_css, J_IDs, J_vals = GISSMOReader.loadcouplinginfojson(load_path)

dict_H_IDs_to_g = Dict(H_IDs .=> collect(1:length(H_IDs)))
J_inds = collect( (dict_H_IDs_to_g[J_IDs[i][1]], dict_H_IDs_to_g[J_IDs[i][2]]) for i = 1:length(J_IDs) )
H_inds = collect( dict_H_IDs_to_g[H_IDs[i]] for i = 1:length(H_IDs))




@assert 1==2

# others.
Id, Ix, Iy_no_im, Iz, Ip, Im, Iy,
    Im_full, Iz_full, Ix_full, Iys_no_im_full = NMRSpectraSimulator.prepcouplingalgorithm(3)


s, Q = eigen(A+B)

sr, Qr = eigen(A)

import Distances
#dist = Distances.Hamming()
dist = Distances.Euclidean()
d = Distances.colwise(dist, Q, Qr)
d2 = Distances.pairwise(dist, Q, Qr, dims=2)

@assert 1==2

tol_coherence = 1e-2

Id, Ix, Iy_no_im, Iz, Ip, Im, Iy,
        Im_full, Iz_full, Ix_full, Iys_no_im_full = NMRSpectraSimulator.prepcouplingalgorithm(length(H_IDs))


# spin Hamiltonian
ω0 = ppm2hzfunc.(H_css)
H, H1, H2 = NMRSpectraSimulator.getgenericHamiltoniantest( Id,
    Ix,
    Iy_no_im,
    Iz,
    ω0,
    J_vals,
    J_inds)

#s_H, Q_H = eigen(H)

a_H, F_H, p_H, s_H, Q_H, coherence_labels_H, M_array_H = NMRSpectraSimulator.getaΩ(Iz_full, H,
        Iys_no_im_full, Ix_full; tol = tol_coherence)

###


# test loading.
load_path = "/home/roy/MEGAsync/inputs/NMR/debug/ethanol_eq.json"
#load_path = "/home/roy/MEGAsync/inputs/NMR/debug/ethanol2.json"
H_IDs, H_css, J_IDs, J_vals = GISSMOReader.loadcouplinginfojson(load_path)

dict_H_IDs_to_g = Dict(H_IDs .=> collect(1:length(H_IDs)))
J_inds = collect( (dict_H_IDs_to_g[J_IDs[i][1]], dict_H_IDs_to_g[J_IDs[i][2]]) for i = 1:length(J_IDs) )
H_inds = collect( dict_H_IDs_to_g[H_IDs[i]] for i = 1:length(H_IDs))

ω0 = ppm2hzfunc.(H_css)
G, G1, G2 = NMRSpectraSimulator.getgenericHamiltoniantest( Id,
    Ix,
    Iy_no_im,
    Iz,
    ω0,
    J_vals,
    J_inds)
#
#s_G, Q_G = eigen(G)

a_G, F_G, p_G, s_G, Q_G, coherence_labels_G, M_array_G = NMRSpectraSimulator.getaΩ(Iz_full, G,
        Iys_no_im_full, Ix_full; tol = tol_coherence)
