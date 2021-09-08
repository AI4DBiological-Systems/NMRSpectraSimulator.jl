
function forcesymmetric(A::Matrix{T})::Matrix{T} where T <: Real
    return (A+A')./2
end

# get total number of elements.
function getnumelements(X::Vector{Vector{Vector{T}}})::Int where T
    return sum( sum( length(X[n][i]) for i = 1:length(X[n]) ) for n = 1:length(X) )
end

# get the number of entries for n, i dimensions.
function collectnumelements(X::Vector{Vector{Vector{T}}})::Vector{Vector{Int}} where T
    return collect( collect( length(X[n][i]) for i = 1:length(X[n]) ) for n = 1:length(X) )
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


function runNLopt!(  opt,
                    p0::Vector{T},
                    obj_func,
                    grad_func,
                    p_lbs,
                    p_ubs;
                    max_iters = 10000,
                    xtol_rel = 1e-12) where T
    #
    @assert length(p0) == length(p_lbs) == length(p_ubs)

    opt.maxeval = max_iters
    opt.lower_bounds = p_lbs
    opt.upper_bounds = p_ubs
    opt.xtol_rel = xtol_rel


    opt.min_objective = (xx, gg)->genericNLoptcostfunc!(xx, gg, obj_func, grad_func)

    # optimize.
    (minf, minx, ret) = NLopt.optimize(opt, p0)

    N_evals = opt.numevals

    return minf, minx, ret, N_evals
end

function genericNLoptcostfunc!(x::Vector{Float64}, df_x, f, df)::Float64

    #
    if length(df_x) > 0
        df_x[:] = df(x)
    end

    return f(x)
end

function visualizemixtureintervals(metabolite_regions,
    display_inds, P)

    flags_set = Vector{BitVector}(undef, length(display_inds))
    flags_i_set = Vector{Vector{BitVector}}(undef, length(display_inds))

    for j = 1:length(display_inds)
        n_select = display_inds[j]
        Q = metabolite_regions[n_select]
        st = Q.interval_starts
        fin = Q.interval_fins

        flags_set[j] = falses(length(P))
        flags_i_set[j] = Vector{BitVector}(undef, length(st))

        for i = 1:length(st)
            flags_i = st[i] .< P .< fin[i]

            # debug.
            flags_i_set[j][i] = flags_i

            flags_set[j] = flags_set[j] .| flags_i
        end
    end

    return flags_set, flags_i_set
end



### from Utilities.jl.

"""
    convertcompactdomain(x::T, a::T, b::T, c::T, d::T)::T

converts compact domain x ∈ [a,b] to compact domain out ∈ [c,d].
"""
function convertcompactdomain(x::T, a::T, b::T, c::T, d::T)::T where T <: Real

    return (x-a)*(d-c)/(b-a)+c
end

function isnumericallyclose(x::T, y::T, tol = eps(T)*2) where T
    if abs(x-y) < tol
        return true
    end

    return false
end
