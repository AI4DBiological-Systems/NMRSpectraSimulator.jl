"""
postprocessmageq(mag_eq_inds::Vector{Vector{Int}}, dict_H_inds_to_css;
    atol = 1e-6)

`mag_eq_inds` contains the magnetic equivalent nuclei indices for a spin system.
This function combine the groups of equivalent nuclei if they have the same chemical shift.

mag_eq_inds[k][l] is the k-th equivalent nuclei, l-th nuclei index.
if mag_eq_inds is not empty, then the array mag_eq_inds[k] must not be empty, for all k.
"""
function postprocessmageq(mag_eq_inds::Vector{Vector{Int}},
    dict_H_inds_to_css;
    atol = 1e-6)

    if isempty(mag_eq_inds)
        return Vector{Vector{Int}}(undef, 0)
    end

    cs = collect( dict_H_inds_to_css[mag_eq_inds[k][1]] for k = 1:length(mag_eq_inds) )

    inds = mag_eq_inds

    # check if all cs are the same.
    if abs(sum(cs) / length(cs) - cs[1]) < atol

        inds, cs_unique = NMRSpectraSimulator.removeredundantsinglets(mag_eq_inds, cs)

        inds = collect( unique(inds[i]) for i = 1:length(inds) )
    end

    return inds
end