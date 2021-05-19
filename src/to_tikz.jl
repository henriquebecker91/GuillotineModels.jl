"""
    to_tikz_picture(p :: CutPattern{D, S}) :: String

Creates a `String` with tikz code representing the `CutPattern`.

This is the same as `to_tikz_rectangles` but wraps the rectangles inside
a "`\\begin{tikzpicture}` ... `\\end{tikzpicture}`" environment
(i.e., the first and last lines are environment declarations).
"""
function to_tikz_picture(p :: CutPattern{D, S}) :: String where {D, S}
	return "\\begin{tikzpicture}\n$(to_tikz_rectangles(p))\n\\end{tikzpicture}"
end

"""
    to_rectangles(p :: CutPattern{D, S}[, x, y])

Returns a `Vector` with a normalized spatial positioning of pieces and plates.

The returned value is a `Vector{Tuple{D, NTuple{4, S}}}`, the first value
is the piece/plate `piece_idx`, the second is a quadruplet `(xo, yo, xf, yf)`
(i.e., the origins and the final coordinates in each dimension).
"""
function to_rectangles(
	p :: CutPattern{D, S}, x :: S = zero(S), y :: S = zero(S)
) where {D, S}
	outer_rectangle = (p.piece_idx, (x, y, x + p.length, y + p.width))

	if !iszero(p.piece_idx)
		@assert isempty(p.subpatterns)
		return [outer_rectangle]
	end

	inner_rectangles = Tuple{D, NTuple{4, S}}[]

	x_, y_ = x, y
	for subpattern in p.subpatterns
		append!(inner_rectangles, to_rectangles(subpattern, x_, y_))
		if p.cuts_are_vertical
			y_ += subpattern.width
		else
			x_ += subpattern.length
		end
	end

	return push!(inner_rectangles, outer_rectangle)
end

"""
    to_tikz_rectangles(p :: CutPattern{D, S}) :: String

Creates a `String` with tikz `\\draw` commands representing the CutPattern.

The tikz code has nodes indicating the piece indexes (waste pieces are not
labeled because we do not distinguish between the waste pieces and the
intermediary plates yet). The same node command is also outputted commented and
with the size of piece, for easy alternance between these two informations.
"""
function to_tikz_rectangles(
	p :: CutPattern{D, S}# future keyword arguments
) :: String where {D, S}
	return join(map(to_rectangles(p)) do rectangle
		(i, (xo, yo, xf, yf)) = rectangle
		s = "\\draw[dashed, thick, black] ($xo, $yo) rectangle ($xf, $yf);"
		if !iszero(i)
			cl, cw = xo + (xf - xo)/2, yo + (yf - yo)/2
			s = s * "\n\\node[font=\\LARGE] at ($cl, $cw) {$i};"
			s = s * "\n%\\node[font=\\LARGE] at ($cl, $cw) {$(xf - xo)x$(yf - yo)};"
		end
		s
	end, '\n')
end

