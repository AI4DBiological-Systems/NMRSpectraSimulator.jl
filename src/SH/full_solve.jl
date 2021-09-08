#### all methods related to fullsolve but non-updated coherence.




# mutates r, A, F.
function setupapproximateΩfullsolve!( r::Vector{T}, # eigenvalue buffer.
                                    A::Vector{T},
                                    F::Vector{T},
                                    states, state_labels,
                        H0,
                        Qs0,
                        fs,
                        SW;
                        update_α_flag::Bool = false) where T <: Real

    #
    N_spins = round(Int, log(2,size(H0,1)))

    Ix = getsingleIx()
    Iy_no_im = getsingleIynoim()
    Ip = getsingleI⁺()
    Iz = getsingleIz()

    Ips = collect( Ip for j = 1:N_spins )
    #Ip_full = Kronecker.kroneckersum(Ips...)
    Ip_full = Kronecker.kron(Ips...)

    Iys_no_im = collect( Iy_no_im for j = 1:N_spins )
    Iys_no_im_full = Kronecker.kroneckersum(Iys_no_im...)

    Ixs = collect( Ix for j = 1:N_spins )
    Ix_full = Kronecker.kroneckersum(Ixs...)

    updateΩfunc = ddcc->getΩfulleigen!(r, F, A, ddcc, H0, fs, SW,
                            Iz, states, Qs0, state_labels,
                            Iys_no_im_full, Ix_full, Ip_full;
                            update_α_flag = update_α_flag)

    return updateΩfunc
end

# mutates r and F and A.
# input is Δcs, which is not Δp_cs.
function getΩfulleigen!(r::Vector{T},
                        F::Vector{T},
                        A,
                        Δcs::Vector{T},
                        H0::Matrix{T},
                        fs::T,
                        SW::T,
                        Iz,
                        states::Vector{Int},
                        Qs0::Vector{Vector{T}},
                        state_labels,
                        Iys_no_im_full,
                        Ix_full,
                        Ip_full;
                        update_α_flag::Bool = false) where T <: Real

    ### set up.
    N_spins = length(Δcs)
    N = length(Qs0)
    resize!(r, N)
    #fill!(r, Inf) # debug.

    ### get eigenpair of new H.
    ΔD = constructΔH(Δcs, fs, SW, Iz)

    # #### debug.
    # cs = Δcs
    # ω0 = (4029.9605179550854 .+ ((cs) .* (fs/SW))) .* 2*π
    # ωIzs = collect( ω0[j] .* Iz for j = 1:length(cs) )
    # ΔD = Kronecker.kroneckersum(ωIzs...)
    # #### end debug.

    H = ΔD + H0 # inefficient.
    #println("norm(ΔD) = ", norm(ΔD))

    #s, Q = eigen(H0) # debug.
    s, Q = eigen(H)


    ### update the eigenvalues r for select states.
    # This is essentially a matching via alignment of eigenvectors to the given states.

    Qs = Vector{Vector{T}}(undef, size(H,1))

    for i = 1:length(states)

        k = states[i]
        target_q = Qs0[k]


        # find the eigenvector of H (which are the columns of Q) that is
        #   closest (in l-2 norm) to the k-th eigenvector of H0
        #   (i.e. Qs0[k]).
        # dists = collect( norm( Q[:,j] - target_q ) for j = 1:size(Q,2) )
        # _, ind = findmin(dists)

        dists = collect( abs(dot( Q[:,j], target_q )) for j = 1:size(Q,2) )
        _, ind = findmax(dists)

        # println("i = ", i)
        # println("k = ", k)
        # println()

        r[k] = s[ind]
        Qs[k] = target_q
    end

    # r[:] = s # debug

    ### get the frequencies.
    L = length(state_labels)
    resize!(F, L)

    for l = 1:L
        F[l] = r[state_labels[l][2]] - r[state_labels[l][1]]
    end

    # sol = findfirst(xx->xx==82, states)
    # println("sol = ", sol)

    if update_α_flag

        # pre-compute.
        IyQs = Vector{Vector{T}}(undef, size(Q,1))
        IxQs = Vector{Vector{T}}(undef, size(Q,1))

        for i in states
            IyQs[i] = Iys_no_im_full*Qs[i] # no im since official formula for α is div by 2im.
            IxQs[i] = Ix_full*Qs[i]
        end

        resize!(A, L)
        for l = 1:L
            ind_s = state_labels[l][2]
            ind_r = state_labels[l][1]

            # println("ind_r = ", ind_r)
            # println("ind_s = ", ind_s)

            t1 = -dot(conj(Qs[ind_r]), IyQs[ind_s])
            t2 = dot(conj(Qs[ind_s]), IxQs[ind_r])

            A[l] = 2*(t1 * t2)

            # # try shift operator.
            # A[l] = dot(conj(Qs[ind_r]), Ip_full*Qs[ind_s])
        end

        normalizetoNspins!(A, N_spins)
    end

    return nothing
    #return s, Q, H# debug.
end

function constructΔH(Δω0::Vector{T}, #Δcs::Vector{T},
                    fs::T,
                    SW::T,
                    Iz) where T

    #Δω0 = Δcs*2*π*fs/SW

    ωIzs = collect( Δω0[j] .* Iz for j = 1:length(Δω0) )
    ΔH = Kronecker.kroneckersum(ωIzs...)

    return ΔH
end
