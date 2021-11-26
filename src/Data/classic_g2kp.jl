struct Classic_G2KP{D, S, P} end

is_collection(_::Val{:Classic_G2KP}) = false
is_collection(_::Classic_G2KP{D, S, P}) where {D, S, P} = false

@auto_hash_equals struct G2KP{D, S, P}
	"Length of the knapsack/'original plate'/'large object'."
	L :: S
	"Width of the knapsack/'original plate'/'large object'."
	W :: S
	"Lengths of the pieces."
	l :: Vector{S}
	"Widths of the pieces."
	w :: Vector{S}
	"Demands for the pieces (i.e., maximum number of copies that may be sold)."
	d :: Vector{D}
	"Profits obtained selling the one unit of the corresponding piece."
	p :: Vector{P}
end

"""
    read_from_string(_ :: Classic_G2KP{D, S, P}, s :: AbstractString)

Expect a string in the following format:
```
<L> <W>
<N>
<l_1> <w_1> <p_1> <d_1>
<l_2> <w_2> <p_2> <d_2>
...
```
`L` and `W` are the original plate length and width. `N` is the number of
pieces in the instance. The `l`, `w`, `p`, and `d`, are the length, width,
profit, and demand, of a piece.

Returns a `G2KP{D, S, P}` object.

# Notes

- If `N` is smaller than the number of lines after it, just the first `N`
  lines will be parsed and passed to the returned object.
- Empty lines and multiples spaces/tabs between columns may happen, but each
  non-empty line should contain the expected number of columns (i.e., the
  first non-empty line should have one number, the second non-empty line two
  numbers, all remaining non-empty lines four numbers).

See also: [`read_from_file`](@ref)
"""
@timeit TIMER function read_from_string(
	_ :: Classic_G2KP{D, S, P}, s :: AbstractString
) :: G2KP{D, S, P} where {D, S, P}
	N = zero(D)
	L = W = zero(S)
	(l, w, p, d) = (
		Vector{S}(), Vector{S}(), Vector{P}(), Vector{D}()
	)
	lnc = one(D)
	for ln in split(s, isequal('\n'); keepempty = false)
		if lnc - 2 == N + 1
			@warn "There are more piece lines than the specified number of pieces." *
				" The extra lines will not be used."
			break
		end
		if      lnc == (1::D)
			(L, W) = map(x->parse(S, x), split(ln))
		elseif  lnc == (2::D)
			N = parse(D, ln)
		else    # all other lines
			piece_cols = split(ln)
			cp_l = parse(S, piece_cols[1])
			cp_w = parse(S, piece_cols[2])
			cp_p = parse(P, piece_cols[3])
			cp_d = parse(D, piece_cols[4])
			push!.((l, w, p, d), (cp_l, cp_w, cp_p, cp_d))
		end
		# The instances of paper "Improved state space relaxation for
		# constrained two-dimensional guillotine cutting problems" have a N
		# different of the real number of lines because the instances with 50
		# and 25 items are the same, only changing the value (they could
		# have removed the unused lines...).
		lnc = lnc + (1::D)
	end
	return G2KP{D, S, P}(L, W, l, w, d, p)
end

function read_from_string(
	_ :: Val{:Classic_G2KP}, s :: AbstractString
)
	return read_from_string(Classic_G2KP{Int, Int, Int}(), s)
end

# === `GuillotineModels.solution_value` implementations

import ..CutPattern
function solution_value(
	:: Val{:G2KP}, instance :: G2KP{D, S, P}, cp :: CutPattern{D, S}
) where {D, S, P}
	p = instance.p
	if iszero(cp.piece_idx)
		z = zero(eltype(p))
		for sp in cp.subpatterns
			z += solution_value(Val(:G2KP), instance, sp)
		end
		return z
	else
		return p[cp.piece_idx]
	end
end

