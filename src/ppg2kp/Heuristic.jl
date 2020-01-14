module Heuristic

using Random

export first_fit_decr, iterated_greedy

# Returns zero if no slice 1:k (for any k, but without changing the ordering)
# exists such that sum(p[1:k]) > bkv and sum(a[1:k]) <= A; otherwise, returns
# the smallest k value which satisfy these conditions.
function promising_kfirst(
	::Type{D},
	p :: AbstractVector{P},
	a :: AbstractVector{P},
	bkv :: P,
	A :: P
) :: D where {D, P}
	@assert length(p) == length(a)
	sump, suma = zero(P), zero(P)
	for k = 1:length(p)
		suma += a[k]
		suma > A && return zero(D)
		sump += p[k]
		sump > bkv && return convert(D, k)
	end
	# the line below is only executed if sum(p) <= bkv
	return zero(D)
end

# Takes two vectors of the same size. Expand the second based on the values
# of the first (i.e., parameter 'd'). Ex.: expand([0, 1, 2, 3], [4, 5, 6, 7])
# gives [5, 6, 6, 7, 7, 7].
function expand(
	d :: AbstractVector{D}, a :: AbstractVector{T}
) :: AbstractVector{T} where {D, T}
	n = length(d)
	@assert length(a) == n
	collect(Base.Iterators.Flatten(
		map(zip(a, d)) do (v, q)
			vs = similar(a, 1)
			vs[1] = v
			repeat(vs, q)
		end
	))
end

# Note: the rng object is modified, as it would be expected.
function iterated_greedy(
	d :: AbstractVector{D},
	p :: AbstractVector{P},
	l :: AbstractVector{S},
	w :: AbstractVector{S},
	L :: S,
	W :: S,
	rng,
	max_iter_since_last_improv :: I = 1000000
) :: Tuple{P, Vector{D}, Vector{Vector{D}}} where {D, S, P, I}
	n = length(d) # number of distinct piece types
	A = convert(P, L) * convert(P, W)
	# assert explanation: all vectors should be the same length.
	@assert [n] == unique!(length.([p, l, w]))
	# loop control variables
	curr_iter = 0 # current iteration
	iter_last_improv = 0 # last iteration in which bkv changed
	# The triple from first_fit_decr (they will be overwritten unless n == 0)
	bkv = zero(P) # best known value
	sel = zeros(D, n) # like 'd' but with the pieces selected for the solution
	pat = Vector{Vector{D}}() # the cut pattern (is always a 2-stage cut)
	# Expand all vector by the demand. This corrects the algorithm for instances
	# in which not all piece types have just a single piece.
	pe, le, we = expand.((d, d, d), (p, l, w))
	# Create a first ordering by efficiency (i.e., cost-benefit, profit-to-area
	# ratio). This way, the same instance with the same RNG will have the same
	# results, independent from the initial ordering of the piece types.
	ie = expand(d, collect(1:n)) # permutation used to keep the original order
	ae = @. convert(P, le) * convert(P, we)
	ee = @. pe / convert(Float64, ae)
	oe = sortperm(ee)

	while curr_iter - iter_last_improv <= max_iter_since_last_improv
		curr_iter += 1
		permute!.((pe, le, we, ae, ie), (oe, oe, oe, oe, oe))
		k = promising_kfirst(D, pe, ae, bkv, A)
		if !iszero(k)
			@views pe_,le_,we_,ae_,ie_ = pe[1:k],le[1:k],we[1:k],ae[1:k],ie[1:k]
			oe_ = sortperm(we_; rev = true)
			permute!.((pe_, le_, we_, ae_, ie_), (oe_, oe_, oe_, oe_, oe_))
			new_v, new_sel, new_pat = first_fit_decr(pe_, le_, we_, ie_, n, L, W)
			@assert sum(new_sel .* p) == new_v
			if new_v > bkv
				bkv, sel, pat = new_v, new_sel, new_pat
			end
		end
		shuffle!(rng, oe)
	end

	nz_i = collect(i for i in 1:n if !iszero(sel[i]))
	@show nz_i
	nz_sel = collect(sel[i] for i in 1:n if !iszero(sel[i]))
	@show nz_sel
	nz_l = collect(l[i] for i in 1:n if !iszero(sel[i]))
	@show nz_l
	nz_w = collect(w[i] for i in 1:n if !iszero(sel[i]))
	@show nz_w
	return bkv, sel, pat
end

function first_fit_decr(
	p :: AbstractVector{P}, # profit of the piece
	l :: AbstractVector{S}, # length of the piece
	w :: AbstractVector{S}, # width of the piece
	t :: AbstractVector{D}, # type of the piece (index at the original demand)
	n :: D, # number of distinct piece types (n >= maximum(t))
	L :: S, # length of the plate
	W :: S # width of the plate
) :: Tuple{P, Vector{D}, Vector{Vector{D}}} where {D, S, P}
	# assert explanation: all vectors should be the same length.
	@assert isone(length(unique!(length.([p, l, w, t]))))
	# assert explanation: this is first-fit WIDTH-decreasing
	@assert issorted(w; rev = true)
	lw = Vector{S}() # levels width
	ll = Vector{S}() # levels length
	rw = W # remaining width

	bkv = zero(P) # best known value
	sel = zeros(D, n) # if piece i is cut then sel[t[i]] is incremented
	pat = Vector{Vector{D}}() # the pattern (always a 2-stage pattern)
	for i in 1:length(p) # for each piece
		assigned = false
		for j in 1:length(lw) # for each already open level
			# assert reasoning: the levels are created with the width of previous
			# pieces, and the pieces are ordered by decreasing weight, therefore
			# no piece may have more width than an already open level.
			@assert lw[j] >= w[i]
			if ll[j] >= l[i] # checks if the level has enough length yet
				ll[j] -= l[i]
				bkv += p[i]
				sel[t[i]] += one(D)
				push!(pat[j], t[i])
				assigned = true
				break
			end
		end
		# If the piece does not fit any already open level...
		if !assigned
			# ... and there is no space to open a new level, then jump to the
			# next piece type.
			rw < w[i] && continue
			# ... but there is yet space to open a new level, create a new level.
			push!(lw, w[i])
			push!(ll, L - l[i])
			rw -= w[i]
			bkv += p[i]
			sel[t[i]] += one(D)
			push!(pat, [t[i]])
		end
	end

	return bkv, sel, pat
end

end # module

