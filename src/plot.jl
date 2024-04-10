const rna_letters = ["A", "C", "G", "U"]
const dna_letters = ["A", "C", "G", "T"]

@userplot LogoPlot

"""
    logoplot(motif::AbstractMatrix; charnames::AbstractVector{AbstractString}, alphabet_coords=ALPHABET_GLYPHS, do_norm=false, ignore_case=false)
"""
@recipe function f(data::LogoPlot; 
    alphabet_coords=ALPHABET_GLYPHS, 
    do_norm=false, ignore_case=false, rna=true)
    motif = data.args[1]
    # charnames = data.args[2]
    charnames= append!(rna ? rna_letters : dna_letters, ["$i" for i in 1:size(motif,1)-4])

    length(charnames) != size(motif, 1) && throw(ArgumentError("number of rows in `motif` matrix does not match length of `charnames`, motif has $(size(motif,1)) rows and charnames has length $(length(charnames))"))
    
    if do_norm
        motif = motif ./ sum(motif, dims=1)
    end

    if ignore_case
        charnames = uppercase.(charnames)
    end

    all_coords = rna ? transform_probs_to_coords_crosslink(motif; charnames, alphabet_coords) : transform_probs_to_coords(motif, charnames; alphabet_coords)
    xticks_pos = size(motif,2) < 5 ? (1:size(motif,2)) : [1, (0:5:size(motif,2))[2:end]...]
    
    grid --> false 
    tickdir --> :out
    # framestyle --> :box
    legend --> false
    size --> logo_size
    xticks --> xticks_pos
    ylabel --> "bits"
    # ylims --> (-crosslink_stretch_factor,2.0)
    thickness_scaling --> 0.0
    # legend --> :outerright
    # margin --> Plots.mm
    for (k, v) in all_coords
        @series begin
            fill := 0
            lw --> 0
            label --> k
            color --> get(AA_PALETTE3, k, :grey)
            v.xs, v.ys
        end
    end
end


