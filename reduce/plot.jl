function batchplot(As, Bs, save_root_folder, molecule_names, compound_select,
    display_reduction_factor, display_threshold_factor;
    canvas_size = (1000, 400),
    display_flag = false)

    A = As[compound_select]
    B = Bs[compound_select]
    name = molecule_names[compound_select]

    # create the functions for each resonance group.
    qs = collect( collect( ωω->B.qs[i][k](ωω-B.ss_params.d[i], B.ss_params.κs_λ[i]) for k = 1:length(B.qs[i]) ) for i = 1:length(B.qs) )
    q_singlets = ωω->NMRSpectraSimulator.evalsinglets(ωω, B.d_singlets, A.αs_singlets, A.Ωs_singlets, B.β_singlets, B.λ0, B.κs_λ_singlets)

    # create the function for the entire compound.
    q = uu->NMRSpectraSimulator.evalitpproxymixture(uu, As[compound_select:compound_select], Bs[compound_select:compound_select])

    # evaluate at the plotting positions.
    q_U = q.(U_rad)

    qs_U = collect( collect( qs[i][k].(U_rad) for k = 1:length(qs[i]) ) for i = 1:length(qs) )
    q_singlets_U = q_singlets.(U_rad)


    # sanity check.
    q_check_U = q_singlets_U
    if !isempty(qs) # some compounds only have singlets.
        q_check_U += sum( sum( qs[i][k].(U_rad) for k = 1:length(qs[i]) ) for i = 1:length(qs) )
    end
    discrepancy = norm(q_check_U- q_U)
    println("sanity check. This should be numerically zero: ", discrepancy)

    # reduce the plotting positions for low signal regions. Otherwise the plot store size will be too large, and the time to load the plot will be long.
    inds, _ = NMRSpectraSimulator.prunelowsignalentries(q_U, display_threshold_factor, display_reduction_factor)
    P_display = P[inds]
    U_display = U[inds]

    # plot.
    save_folder = joinpath(save_root_folder, molecule_names[compound_select])
    isdir(save_folder) || mkpath(save_folder)

    plots_save_path = joinpath(save_folder, "groups_real.html")
    title_string = "$(name) resonance groups, real part"
    plot_obj, q_U, qs_U, q_singlets_U = plotgroups(title_string, P_display, U_display, q, qs, q_singlets, real, P[1]; canvas_size = canvas_size)
    Plots.savefig(plot_obj, plots_save_path)

    if display_flag
        display(plot_obj)
    end

    return nothing
end

"""
save resonance groupings of a compound.
Choices for `f`` are: `real()`, `imag()`, or `abs()`. These corresponds to real part, imaginary part, and magnitude spectrum, respectively.
"""
function plotgroups(title_string::String,
    P,
    U,
    q,
    qs,
    q_singlets,
    f::Function,
    return_val_type::T;
    canvas_size::Tuple{Int,Int} = (1000, 400)) where T

    U_rad = U .* (2*π)

    q_U = q.(U_rad)
    plot_obj = Plots.plot( P,
        f.(q_U),
        title = title_string,
        label = "sum of sub-models",
        seriestype = :line,
        ticks = :native,
        xlims = (P[1],P[end]),
        hover = P,
        linewidth = 4,
        xflip = true,
        size = canvas_size)

    qs_U = Vector{Vector{Vector{Complex{T}}}}(undef, length(qs))
    for i = 1:length(qs)

        qs_U[i] = Vector{Vector{Complex{T}}}(undef, length(qs[i]))
        for k = 1:length(qs[i])

            qs_U[i][k] = qs[i][k].(U_rad)

            Plots.plot!(plot_obj, P, f.(qs_U[i][k]), label = "sub-group ($(i),$(k))",
                seriestype = :line,
                linestyle = :dot,
                xflip = true,
                linewidth = 4)
        end
    end

    q_singlets_U = q_singlets.(U_rad)
    Plots.plot!(plot_obj, P, f.(q_singlets_U), label = "group of all singlets",
        seriestype = :line,
        linestyle = :dot,
        xflip = true,
        linewidth = 4)
    #
    #Plots.plot!(plot_obj, P, f.(q_U),
    #markershape = :circle,
    #seriestype = :scatter,
    #xflip = true)

    return plot_obj, q_U, qs_U, q_singlets_U
end





println("Timing fitproxies()")
@time Bs = NMRSpectraSimulator.fitproxies(As, dummy_SSFID, λ0;
    #names = molecule_names[compound_select:compound_select],
    names = molecule_names,
    config_path = surrogate_config_path,
    Δcs_max_scalar_default = Δcs_max_scalar_default,
    κ_λ_lb_default = κ_λ_lb_default,
    κ_λ_ub_default = κ_λ_ub_default,
    u_min = u_min,
    u_max = u_max,
    Δr_default = Δr_default,
    Δκ_λ_default = Δκ_λ_default)

for compound_select = 1:length(As)

    display_reduction_factor = 100
    display_threshold_factor =  0.05/10

    batchplot(As, Bs, save_root_folder, molecule_names, compound_select,
        display_reduction_factor, display_threshold_factor;
        canvas_size = (1000, 400),
        display_flag = false)
end

# display number of resonance groups per spin system that aren't singlets.
for n = 1:length(As)
    A = As[n]
    N_β_vars_sys = collect( length(A.Δc_bar[i][1]) for i = 1:length(A.Δc_bar) )

    println("$(n). $(molecule_names[n]) N_β_vars_sys: ", N_β_vars_sys)
end
