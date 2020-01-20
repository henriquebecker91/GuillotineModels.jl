module InstanceReader

using TimerOutputs

export read_from_string, read_from_file

function read_from_string(s :: AbstractString)
	@timeit "read_from_string" begin
	N = L = W = (0 :: Int64)
	(l, w, p, q) = (
		Vector{Int64}(), Vector{Int64}(), Vector{Int64}(), Vector{Int64}()
	)
	lnc = 1
	for ln in split(s, isequal('\n'); keepempty = false)
		if      lnc == 1
			(L, W) = map(x->parse(Int64, x), split(ln))
		elseif  lnc == 2
			N = parse(Int64, ln)
		else    # all other lines
			# The instances of paper "Improved state space relaxation for
			# constrained two-dimensional guillotine cutting problems" have a N
			# different of the real number of lines because the instances with 50
			# and 25 items are the same, only changing the value (they could
			# have removed the unused lines...).
			if (lnc > N + 2)
			end
			# cp_ == current plate
			(cp_l, cp_w, cp_p, cp_q) = map(x -> parse(Int64, x), split(ln))
			push!((l, w, p, q), (cp_l, cp_w, cp_p, cp_q))
		end
		if lnc - 2 == N # use N and not all the lines in the file
			@warn "There are more plate lines than the specified number of plates."
			break
		end
		lnc = lnc + 1
	end
	end # @timeit
	return N, L, W, l, w, p, q, 
end

# TODO: restore the configuration Dict and allow:
# 	to parse instances with an id column; out of order plate ids; save ids
# TODO: SOME INSTANCES ASSUME THE READER IS COMPLETELY INSENSTIVE TO
# 	WHITESPACE, CHANGE THE METHOD TO GET THE NEXT TOKEN AND NOT WORK
# 	ORIENTED BY LINES BUT BY WORDS/TOKENS/NUMBERS
function read_from_file(filepath :: AbstractString)
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

