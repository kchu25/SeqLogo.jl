module SeqLogo

using Plots
# Write your package code here.
include("glyphs.jl")
include("transform.jl")
include("plot.jl")

logo_size = (600, 200)
letter_scale = 2
crosslink_stretch_factor = 9

function save_crosslinked_pwm(motif::AbstractMatrix, save_where::String)
    @assert size(motif,1) ≥ 5 "crosslinked motif must have at least 5 rows"
    p = logoplot(motif; do_norm=false, ignore_case=false);
    savefig(p, save_where)
end


end
