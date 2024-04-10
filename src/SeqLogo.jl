module SeqLogo

using Plots
# Write your package code here.
include("glyphs.jl")
include("transform.jl")
include("plot.jl")

logo_size = (800, 200)
letter_scale = 2
crosslink_stretch_factor = 9

function save_crosslinked_pwm(motif::AbstractMatrix, save_where::Union{String, Nothing} = nothing)
    @assert size(motif,1) ≥ 5 "crosslinked motif must have at least 5 rows"
    p = logoplot(motif; do_norm=false, ignore_case=false);
    if isnothing(save_where)
        return p
    end
    savefig(p, save_where)
end

function plot_pwm(motif::AbstractMatrix, save_where::Union{String, Nothing} = nothing)
    @assert size(motif,1) ≥ 4 "motif must have at least 4 rows"
    p = logoplot(motif; do_norm=false, ignore_case=false, rna=false);
    if isnothing(save_where)
        return p
    end
    savefig(p, save_where)
end



end
