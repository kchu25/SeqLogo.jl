module SeqLogo

using Plots
# Write your package code here.
include("glyphs.jl")
include("transform.jl")
include("plot.jl")

crosslink_stretch_factor = 10

function save_crosslinked_pwm(motif::AbstractMatrix, save_where::String)
    @assert size(motif,1) == 5 "crosslinked motif must have 5 rows"
    p = logoplot(motif; do_norm=false, ignore_case=false);
    savefig(p, save_where)
end


end
