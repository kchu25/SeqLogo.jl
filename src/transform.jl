
ic_height_uniform(x; bg=0.25, ϵ=1e-20) = (x+ϵ) * log2( (x+ϵ) / bg)
ic_height_here(col) = sum(ic_height_uniform.(col))

function transform_probs_to_coords(motif, charnames; alphabet_coords=ALPHABET_GLYPHS)
    all_coords = []
    # for each character (row) collect all positions and heights of that character's polygon
    for (j, c) in enumerate(charnames)
        xs, ys = Float64[], Float64[]
        # get character glyph coords else fall back to a simple rectangle
        charglyph = get(alphabet_coords, c, BASIC_RECT)
        # for each postion in the sequence push in the coords for the character's polygon and adjust y_height based on the motif's weight
        for  (xoffset,col) in enumerate(eachcol(motif))
            ic_height = ic_height_here(col)
            adjusted_height = ic_height .* col
            yoffset = sum(adjusted_height[adjusted_height .< adjusted_height[j]])
            push!(xs, (1.2 .* charglyph.x .+ xoffset .- .5)...)
            push!(xs, NaN)
            push!(ys, (adjusted_height[j] .* charglyph.y .+ yoffset)...)
            push!(ys, NaN)
        end
        push!(all_coords, (c, (;c=[c], xs, ys)))
    end
    all_coords
end

function transform_probs_to_coords_crosslink(motif; charnames=["A","C","G","T","*"], alphabet_coords=ALPHABET_GLYPHS)
    # "" means returning BASIC_RECT
    all_coords = []
    # for each character (row) collect all positions and heights of that character's polygon
    for (j, c) in enumerate(charnames)
        xs, ys = Float64[], Float64[]
        # get character glyph coords else fall back to a simple rectangle
        charglyph = get(alphabet_coords, c, BASIC_RECT)
        # for each postion in the sequence push in the coords for the character's polygon and adjust y_height based on the motif's weight
        if charglyph == BASIC_RECT
            for (xoffset,col) in enumerate(eachcol(motif))
                push!(xs, (1.2 .* charglyph.x .+ xoffset .- .5)...)
                push!(xs, NaN)
                # yoffset_c = 0.0
                adjusted_height = col[5:end] .* crosslink_stretch_factor
                total_height = sum(adjusted_height)
                yoffset = sum(adjusted_height[adjusted_height .< adjusted_height[j-4]])
                push!(ys, (adjusted_height[j-4] .* charglyph.y .- (total_height - yoffset))...)
                push!(ys, NaN)
            end
        else
            for (xoffset,col) in enumerate(eachcol(motif))
                acgt = @view col[1:4]
                ic_height = ic_height_here(acgt)
                adjusted_height = ic_height .* acgt
                # println("$(adjusted_height .< adjusted_height[j])")
                # println("$(sum(adjusted_height[adjusted_height .< adjusted_height[j]]))")
                yoffset = sum(adjusted_height[adjusted_height .< adjusted_height[j]])
                push!(xs, (1.2 .* charglyph.x .+ xoffset .- .5)...)
                push!(xs, NaN)
                push!(ys, letter_scale .* (adjusted_height[j] .* charglyph.y .+ yoffset)...)
                push!(ys, NaN)
            end
        end
        push!(all_coords, (c, (;c=[c], xs, ys)))
    end
    all_coords
end
