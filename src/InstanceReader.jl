"""
InstanceReader groups methods that parse instance files.
For now, only the knapsack demand-constrained variant format is supported, and
only without the id column (each plate should have only four values: length,
width, profit, and demand).
"""
module InstanceReader

using TimerOutputs

export read_from_string, read_from_file

"""
    read_from_string(s :: AbstractString)

Expect a string in the following format:
```
<L> <W>
<N>
<l_1> <w_1> <p_1> <d_1>
<l_2> <w_2> <p_2> <d_2>
...
```
Returns a tuple of `(N, L, W, l, w, p, d) :: Tuple{Int64, Int64, Int64,
Vector{Int64}, Vector{Int64}, Vector{Int64}, Vector{Int64}}`.
`L` and `W` are the original plate length and width. `N` is the number of
pieces in the instance. The `l`, `w`, `p`, and `d`, are the length, width,
profit, and demand, of a piece.

# Notes
- If `N` is smaller than the number of lines after it, just the first `N` lines
will be parsed and added to return tuple.
- For now all values are parsed as Int64, downsizing to a smaller type may
be done after.
- Empty lines and multiples spaces/tabs between columns may happen, but each non-empty line should contain the expected number of columns (i.e., the first non-empty line should have one number, the second non-empty line two numbers, all remaining non-empty lines four numbers).

See also: [`read_from_file`](@ref)
"""
function read_from_string(
	s :: AbstractString
) :: Tuple{
	Int64,Int64,Int64,Vector{Int64},Vector{Int64},Vector{Int64},Vector{Int64}
}
	@timeit "read_from_string" begin
	N = L = W = (0 :: Int64)
	(l, w, p, q) = (
		Vector{Int64}(), Vector{Int64}(), Vector{Int64}(), Vector{Int64}()
	)
	lnc = 1
	for ln in split(s, isequal('\n'); keepempty = false)
		if lnc - 2 == N + 1 # use N and not all the lines in the file
			@warn "There are more plate lines than the specified number of plates." *
				" The extra lines will not be used."
			break
		end
		if      lnc == 1
			(L, W) = map(x->parse(Int64, x), split(ln))
		elseif  lnc == 2
			N = parse(Int64, ln)
		else    # all other lines
			# cp_ == current plate
			(cp_l, cp_w, cp_p, cp_q) = map(x -> parse(Int64, x), split(ln))
			push!.((l, w, p, q), (cp_l, cp_w, cp_p, cp_q))
		end
		# The instances of paper "Improved state space relaxation for
		# constrained two-dimensional guillotine cutting problems" have a N
		# different of the real number of lines because the instances with 50
		# and 25 items are the same, only changing the value (they could
		# have removed the unused lines...).
		lnc = lnc + 1
	end
	end # @timeit
	return N, L, W, l, w, p, q
end

# TODO: restore the configuration Dict and allow:
# 	to parse instances with an id column; out of order plate ids; save ids
# TODO: SOME INSTANCES ASSUME THE READER IS COMPLETELY INSENSTIVE TO
# 	WHITESPACE, CHANGE THE METHOD TO GET THE NEXT TOKEN AND NOT WORK
# 	ORIENTED BY LINES BUT BY WORDS/TOKENS/NUMBERS
"""
    read_from_file(filepath :: AbstractString)

Reads the given file and apply `read_from_string` to its contents.

See also: [`read_from_string`](@ref)
"""
function read_from_file(
	filepath :: AbstractString
) :: Tuple{
	Int64,Int64,Int64,Vector{Int64},Vector{Int64},Vector{Int64},Vector{Int64}
}
	@timeit "read_from_file" begin
	# Just declare and define the type, does not matter that is the same shared
	# array, they will refer to others after.
	N = L = W = (0 :: Int64)
	l = w = p = q = Vector{Int64}()
	open(filepath) do f
		N, L, W, l, w, p, q = read_from_string(read(f, String))
	end
	end # @timeit

	return N, L, W, l, w, p, q
end

end # module

