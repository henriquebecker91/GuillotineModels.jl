"""
GuillotineModels.PPG2KP.Heuristic implements the heuristic needed for faithful
reimplementation of the Priced PP-G2KP from ````F. Furini, E. Malaguti, and D.
Thomopulos, ``Modeling Two-Dimensional Guillotine Cutting Problems via Integer
Programming,'' INFORMS Journal on Computing, vol. 28, no. 4, pp. 736–751, Oct.
2016, doi: 10.1287/ijoc.2016.0710.````. The heuristic itself, however, is older
and described in: ````M. Dolatabadi, A. Lodi, and M. Monaci, ``Exact algorithms
for the two-dimensional guillotine knapsack,'' Computers \\& Operations
Research, vol. 39, no. 1, pp. 48–53, Jan. 2012, doi: 10.1016/j.cor.2010.12.018.
````.

The methods exported by this module are not of general interest, you need to
either want: (1) a fast but shelf-restricted heuristic for the guillotine 2D
knapsack problem; (2) to check if they were implemented correctly;
(3) to use them to make your own implementation of Priced PP-G2KP.
"""
module Heuristic

using Random

export first_fit_decr, iterated_greedy

"""
    promising_kfirst(::Type{D}, p, a, bkv, A) :: D

Zero if there is no `k` for which `sum(p[1:k]) > bkv && sum(a[1:k]) <= A`;
otherwise, return the smallest `k` value satisfying these conditions.

In other words, without changing the order of the given sequence, return size
of the prefix of smallest size (a piece subset) that has better profit than the
best known value and that it is not immediately obvious that it cannot fit the
given area (does not break the naive area bound). If such prefix does not
exist, then just return zero.
"""
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

"""
    expand(d :: [D], a :: [T]) :: [T]

Given two vectors of the same size, create a copy of `a` that replaces each
element `a[i]` by `d[i]` copies of it, and then flatten the copy.

```julia
expand([0, 1, 2, 3], [4, 5, 6, 7])
[5, 6, 6, 7, 7, 7]
```
"""
function expand(
	d :: AbstractVector{D}, a :: AbstractVector{T}
) :: Vector{T} where {D, T}
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

"""
    iterated_greedy(d, p, l, w, L, W, rng, max_iter_since_last_improv = 10^6)

The iterated greedy procedure is a heuristic that creates a random piece subset
that may improve the best known value (BKV) and call the first fit
with-decreasing heuristic over it (if this inner heuristic is able to position
all pieces then the BKV is improved).

The first fit heuristic is deterministic and gives the same solution for a
piece subset (the ordering is irrelevant), so the iterated procedure is adding
diversity to it. The first fit heuristic is also shelf-restricted, so the
solution is always shelf-restricted and may be impossible to reach the
optimality (in the case no optimal solution is shelf-restricted).

The procedure is repeated until `max_iter_since_last_improv` calls to
the first fit width-decreasing were made without improving the BKV.

# Arguments

1. `d::AbstractVector{D}`: The pieces demand.
2. `p::AbstractVector{P}`: The pieces profit.
3. `l::AbstractVector{S}`: The pieces length.
4. `w::AbstractVector{S}`: The pieces width.
5. `L::S`: The original plate length.
6. `W::S`: The original plate width.
7. `rng::AbstractRNG`: The random number generator.
8. `max_iter_since_last_improv::I`: The parameter for the stopping criteria,
  after `max_iter_since_last_improv` iterations since the last change in the
  best known value the algorithm ends.

# Returns

The same of the [`first_fit_decr`](@ref) method (the iterated greedy just
returns the best solution found by it).
"""
function iterated_greedy(
	d :: AbstractVector{D},
	p :: AbstractVector{P},
	l :: AbstractVector{S},
	w :: AbstractVector{S},
	L :: S,
	W :: S,
	rng :: AbstractRNG,
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

"""
    first_fit_decr(p, l, w, t, n, L, W) :: (P, [D], [[D]])

The first fit width-decreasing heuristic is a deterministic, fast, and
shelf-restricted heuristic for the knapsack 2D guillotine problem.

This method is mainly called by [`iterated_greedy`](@ref) and because
of this the arguments are adjusted so the method is aware it is receiving a
subset of the complete piece set.

# Arguments

1. `p::AbstractVector{P}`: The pieces profit.
2. `l::AbstractVector{S}`: The pieces length.
3. `w::AbstractVector{S}`: The pieces width.
4. `t::AbstractVector{D}`: The piecex type/index/identifier (index in the
  original piece set).
5. `n::D`: The number of pieces in set (necessary to build the second
  returned vector).
6. `L::S`: The original plate length.
7. `W::S`: The original plate width.

# Returns

1. The best know value, i.e., the value of the returned solution.
2. A vector `v` of the same size as the number of pieces; if `v[i] == x`
  then the piece `i` has `x` copies in the solution returned.
3. A solution with the best known value, as the heuristic only consider
  shelf-restricted configurations the solution may be represented as a
  vector of vectors. Each vector is a length strip (with the width of
	the first piece inside, and the summed length of all pieces inside).
	Inside each vector the pieces are ordered by non-increasing width,
	and the vectors/stripes are sorted by non-increasing width too (
	i.e., the width of their first element).
"""
function first_fit_decr(
	p :: AbstractVector{P},
	l :: AbstractVector{S},
	w :: AbstractVector{S},
	t :: AbstractVector{D},
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
			# pieces, and the pieces are ordered by decreasing width, therefore
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

