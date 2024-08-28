function crosslink!(crosslink_mat, all_coords, logo_x_offset, logo_y_offset)
    if !isnothing(crosslink_mat)
        charglyph = C_RECT;
        num_rows = size(crosslink_mat,1)
        for (xoffset, col) in enumerate(eachcol(crosslink_mat))            
            vsp = 1e-8 .* rand(length(col)) # vsp: very small perturb
            esh = col .* crosslink_stretch_factor .+ vsp # esh: each signiture height
            total_height = sum(esh)
            for r = 1:num_rows
                xs, ys = Float64[], Float64[]
                yoffset = sum(esh[esh .< esh[r]])
                # @info "esh1: $(esh[1]), esh2: $(esh[2])"
                push!(ys, (esh[r] .* charglyph.y .- (total_height - yoffset) .+ logo_y_offset)...)
                push!(ys, NaN)
                push!(xs, ((1.2 .* charglyph.x .+ xoffset .- 0.275)...))
                push!(xs, NaN)
                push!(all_coords, ("$r", (;xs, ys)))
            end
        end
    end
end

function freq2xyWcrosslink(pfm, crosslink_mat; 
                 background=[0.25 for _ = 1:4],
                 rna=true, 
                 beta=1.0, # width of the x-axis for each letter in the motif
                 logo_x_offset = 0.0,
                 logo_y_offset = 0.0, 
                 alphabet_coords=ALPHABET_GLYPHS,
                 very_small_perturb = 1e-5 .* rand(4))
    # @assert sum(pfm, dims=1) .≈ 1 "pfm must be a probability matrix"
    all_coords = []
    charnames = rna ? rna_letters : dna_letters
    # for each character (row) collect all positions and heights of that character's polygon
    for (j, c) in enumerate(charnames)
        xs, ys = Float64[], Float64[]
        # get character glyph coords else fall back to a simple rectangle
        charglyph = get(alphabet_coords, c, BASIC_RECT)
        #= for each postion in the sequence push in the coords 
           for the character's polygon and adjust y_height based 
           on the motif's weight =#
        for (xoffset,col) in enumerate(eachcol(pfm))
            acgt = @view col[1:4]
            ic_height = ic_height_here(col; background=background)
            adjusted_height = ic_height .* acgt .+ very_small_perturb
            yoffset = sum(adjusted_height[adjusted_height .< adjusted_height[j]])
            push!(xs, ((beta * 1.2) .* charglyph.x .+ (1/((beta * 0.9)))*0.35 .+ xoffset .+ (logo_x_offset - 1))...)
            push!(xs, NaN)
            push!(ys, (adjusted_height[j] .* charglyph.y .+ yoffset .+ logo_y_offset)...)
            push!(ys, NaN)
        end        
        push!(all_coords, (c, (;xs, ys)))
    end
    crosslink!(crosslink_mat, all_coords, logo_x_offset, logo_y_offset)
    all_coords
end


@userplot LogoPlotWithCrosslink
@recipe function f(data::LogoPlotWithCrosslink; 
                   rna=true, 
                   xaxis=false,
                   yaxis=false,
                   left_margin=155Plots.mm,
                   margin = 225Plots.mm,
                   ytickfontsize=185,
                   thickness_scaling=0.05,
                   logo_x_offset=0.0,
                   logo_y_offset=0.0,
                   beta=1.0)

    num_cols = size(data.args[1], 2)
    logo_size = (_width_factor_(num_cols)*num_cols, 220)
    framestyle --> :zerolines
    ylims --> (-crosslink_stretch_factor, 2)
    yticks --> 0:1:2  # Ensure ticks are generated
    ytickfontcolor --> :gray
    xtickfontcolor --> :gray
    xticks --> 1:1:num_cols
    # margin --> margin
    left_margin --> left_margin
    ytickfontsize --> ytickfontsize
    xtickfontsize --> 145    
    ytickfont --> font(36, "Helvetica")
    xaxis && (xaxis --> xaxis)
    yaxis && (yaxis --> yaxis)
    legend --> false
    tickdir --> :in
    grid --> false
    thickness_scaling --> thickness_scaling
    size --> logo_size
    pfm = data.args[1]
    background = data.args[2]
    crosslink_mat = data.args[3]
    coords = freq2xyWcrosslink(pfm, crosslink_mat;
                     background=background, rna=rna, beta=beta,
                     logo_x_offset=logo_x_offset,
                     logo_y_offset=logo_y_offset);
    for (k, v) in coords
        @series begin
            fill := 0
            lw --> 0
            label --> k
            color --> get(AA_PALETTE3, k, :grey)
            v.xs, v.ys
        end
    end
end

"""
    save_crosslinked_logoplot(pfm, background, c, save_name; dpi=65, rna=true)
    Save a logoplot with crosslinks to a file
"""
function save_crosslinked_logoplot(pfm, background, c, save_name; dpi=65, rna=true)
    @assert all(sum(pfm, dims=1) .≈ 1) "pfm must be a probability matrix"
    @assert length(background) == 4 "background must be a vector of length 4"
    @assert all(0 .≤ background .≤ 1) "background must be a vector of probabilities"
    @assert sum(background) ≈ 1 "background must sum to 1"
    @assert length(c) == size(pfm, 2) "C must be a vector of length equal to the number of columns in pfm"
    @assert all(0 .≤ c .≤ 1) "C must be a vector of probabilities"
    @assert sum(c) ≤ 1 "The sum of C must be less than or equal 1"
    p = logoplotwithcrosslink(pfm, background, c; dpi=dpi, rna=rna)
    savefig(p, save_name)
end

function save_crosslinked_logoplot(pfm, c, save_name; dpi=65, rna=true)
    p = logoplotwithcrosslink(pfm, [0.25 for _ = 1:4], c; dpi=dpi, rna=rna)
    savefig(p, save_name)
end
