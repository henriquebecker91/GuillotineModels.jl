"""
GuillotineModels.PPG2KP.Heuristic implements the heuristic needed for faithful
reimplementation of the Priced PP-G2KP from \\`F. Furini, E. Malaguti, and D.
Thomopulos, \\`\\`Modeling Two-Dimensional Guillotine Cutting Problems via Integer
Programming,'' INFORMS Journal on Computing, vol. 28, no. 4, pp. 736–751, Oct.
2016, doi: 10.1287/ijoc.2016.0710.' The heuristic itself, however, is older
and described in: \\`M. Dolatabadi, A. Lodi, and M. Monaci, \\`\\`Exact algorithms
for the two-dimensional guillotine knapsack,'' Computers \\& Operations
Research, vol. 39, no. 1, pp. 48–53, Jan. 2012, doi: 10.1016/j.cor.2010.12.018.'

The methods exported by this module are not of general interest, you need to
either want: (1) a fast but shelf-restricted heuristic for the guillotine 2D
knapsack problem; (2) to check if they were implemented correctly;
(3) to use them to make your own implementation of Priced PP-G2KP.
"""
module Heuristic

using Random#: shuffle!
import ...CutPattern
import UnsafeArrays

export first_fit_decr, iterated_greedy
export Levels, fast_first_fit_decr, fast_iterated_greedy

"""
    promising_kfirst(::Type{D}, p, a, bkv, A) :: D

!!! **Internal use.**

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
	for k in eachindex(p)
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

!!! **Internal use.**

Given two vectors of the same size, create a copy of `a` that replaces each
element `a[i]` by `d[i]` copies of it, and then flatten the copy.

```julia
expand([0, 1, 2, 3], [4, 5, 6, 7])
[5, 6, 6, 7, 7, 7]
```
"""
function expand(
	d :: AbstractVector{D}, a :: AbstractVector{T}
) where {D, T}
	n = length(d)
	@assert axes(d) == axes(a)
	sum_d = sum(d)
	b = Vector{T}(undef, sum_d)
	next = 1
	for k in eachindex(d)
		a_k = a[k]
		for _ = 1:d[k]
			b[next] = a_k
			next += 1
		end
	end
	return b
end

"""
    shelves2cutpattern(shelves :: [[D]], l, w, L, W) :: CutPattern{D,S}

Transform the third return of `iterated_greedy` or `first_fit_decr` into a
CutPattern object.
"""
function shelves2cutpattern(
	shelves :: Vector{Vector{D}}, l :: Vector{S}, w :: Vector{S}, L :: S, W :: S
) :: CutPattern{D, S} where {D, S}
	outer_cp = Vector{CutPattern{D, S}}()
	for l_strip in shelves
		inner_cp = Vector{CutPattern{D, S}}()
		for pii in l_strip
			push!(inner_cp, CutPattern(l[pii], w[pii], pii))
		end
		push!(outer_cp, CutPattern(L, w[first(l_strip)], false, inner_cp))
	end
	return CutPattern(L, W, true, outer_cp)
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

	#=
	nz_i = collect(i for i in 1:n if !iszero(sel[i]))
	@show nz_i
	nz_sel = collect(sel[i] for i in 1:n if !iszero(sel[i]))
	@show nz_sel
	nz_l = collect(l[i] for i in 1:n if !iszero(sel[i]))
	@show nz_l
	nz_w = collect(w[i] for i in 1:n if !iszero(sel[i]))
	@show nz_w
	=#
	return bkv, sel, pat
end

# Pre-allocated struct for fast_first_fit_decr.
struct Levels{D, S}
	lw :: Vector{S} # width of the level
	ll :: Vector{S} # length of the level
	la :: Vector{D} # the amount of pieces in the level
	qt :: Ref{Int} # quantity of levels
	function Levels{D, S}(max_length :: Int) where {D, S}
		return new(
			Vector{S}(undef, max_length),
			Vector{S}(undef, max_length),
			Vector{D}(undef, max_length),
			Ref(0)
		)
	end
end

Base.empty!(c :: Levels{D, S}) where {D, S} = c.qt[] = 0

function _extract_shelves(
	pat :: AbstractArray{D, 2},
	lvls :: Levels{D, S}
) where {D, S}
	shelves = [Vector{D}(undef, lvls.la[i]) for i = 1:lvls.qt[]]
	for i = 1:lvls.qt[]
		shelve = shelves[i]
		for j = 1:lvls.la[i]
			shelve[j] = pat[j, i]
		end
	end
	return shelves
end

# CAUTION: this method does not assert it, but the three vectors need to have
# the same axes (not just the length) and all values inside indexes must be
# valid indexes for the two other vectors.
function _swap_permute!(buffer, values, indexes)
	for (new, old) in pairs(indexes)
		buffer[new] = values[old]
	end
	return buffer, values
end

function fast_iterated_greedy(
	d :: AbstractVector{D},
	p :: AbstractVector{P},
	l :: AbstractVector{S},
	w :: AbstractVector{S},
	L :: S,
	W :: S,
	rng :: AbstractRNG,
	max_iter_since_last_improv :: I = 1000000
) :: Tuple{P, Vector{D}, Vector{Vector{D}}} where {D, S, P, I}
	n = convert(D, length(d)) # number of distinct piece types
	A = convert(P, L) * convert(P, W)
	# assert explanation: all vectors should be the same length.
	@assert [n] == unique!(length.([p, l, w]))
	# loop control variables
	curr_iter = 0 # current iteration
	iter_last_improv = 0 # last iteration in which bkv changed
	# Expand all vector by the demand. This corrects the algorithm for instances
	# in which not all piece types have just a single piece.
	pe = expand(d, p)
	le = expand(d, l)
	we = expand(d, w)
	m = length(pe) # the same as sum(d), the number of pieces (not piece types)
	# The values for first_fit_decr (they will be overwritten unless n == 0)
	bkv = zero(P) # best known value
	sel = zeros(D, n) # like 'd' but with the pieces selected for the solution
	pat = Array{D, 2}(undef, (m, m)) # the cut pattern (is always a 2-stage cut)
	lvls = Levels{D, S}(m)
	# Create a first ordering by efficiency (i.e., cost-benefit, profit-to-area
	# ratio). This way, the same instance with the same RNG will have the same
	# results, independent from the initial ordering of the piece types.
	ie = expand(d, collect(one(D):n)) # permutation used to keep the original order
	ae = @. convert(P, le) * convert(P, we)
	ee = @. pe / convert(Float64, ae)
	oe = convert.(D, sortperm(ee))

	tmp_sel = zeros(D, n)
	tmp_pat = Array{D, 2}(undef, (m, m))
	tmp_lvls = Levels{D, S}(m)
	aux_vec_D = Vector{D}(undef, m)
	aux_oe_  = Vector{D}(undef, m)
	aux_vec_S = Vector{S}(undef, m)
	aux_vec_P = Vector{P}(undef, m)
	# I think I do not need UnsafeArrays. The reason for views are:
	# (1) fast_first_fit_decr, that is trivial to receive a `k` value and then
	# just work over that subset; (2) getting the sortperm of a prefix of one
	# of the arrays; (3) permuting the prefix of other arrays inplace considering
	# this prefix.
	UnsafeArrays.@uviews pe le we ae ie aux_oe_ aux_vec_D aux_vec_S aux_vec_P begin
		while curr_iter - iter_last_improv <= max_iter_since_last_improv
			curr_iter += 1
			pe, aux_vec_P = _swap_permute!(aux_vec_P, pe, oe)
			le, aux_vec_S = _swap_permute!(aux_vec_S, le, oe)
			we, aux_vec_S = _swap_permute!(aux_vec_S, we, oe)
			ae, aux_vec_P = _swap_permute!(aux_vec_P, ae, oe)
			ie, aux_vec_D = _swap_permute!(aux_vec_D, ie, oe)
			k = promising_kfirst(D, pe, ae, bkv, A)
			if !iszero(k)
				r = 1:k
				@views pe_, le_, we_, ae_, ie_ = pe[r], le[r], we[r], ae[r], ie[r]
				@views aux_vec_D_, aux_vec_S_, aux_vec_P_ = aux_vec_D[r], aux_vec_S[r], aux_vec_P[r]
				oe_ = @view aux_oe_[r]
				sortperm!(oe_, we_; rev = true)
				permute!.((pe_, le_, we_, ae_, ie_), (oe_, oe_, oe_, oe_, oe_))
				#pe_, aux_vec_P_ = _swap_permute!(aux_vec_P_, pe_, oe_)
				#le_, aux_vec_S_ = _swap_permute!(aux_vec_S_, le_, oe_)
				#we_, aux_vec_S_ = _swap_permute!(aux_vec_S_, we_, oe_)
				#ae_, aux_vec_P_ = _swap_permute!(aux_vec_P_, ae_, oe_)
				#ie_, aux_vec_D_ = _swap_permute!(aux_vec_D_, ie_, oe_)
				new_v = fast_first_fit_decr!(
					tmp_sel, tmp_pat, tmp_lvls, pe_, le_, we_, ie_, n, L, W
				)
				#@assert sum(tmp_sel .* p) == new_v
				if new_v > bkv
					bkv = new_v
					sel, tmp_sel = tmp_sel, sel
					pat, tmp_pat = tmp_pat, pat
					lvls, tmp_lvls = tmp_lvls, lvls
				end
			end
			shuffle!(rng, oe)
		end
	end

	#=
	nz_i = collect(i for i in 1:n if !iszero(sel[i]))
	@show nz_i
	nz_sel = collect(sel[i] for i in 1:n if !iszero(sel[i]))
	@show nz_sel
	nz_l = collect(l[i] for i in 1:n if !iszero(sel[i]))
	@show nz_l
	nz_w = collect(w[i] for i in 1:n if !iszero(sel[i]))
	@show nz_w
	=#
	return bkv, sel, _extract_shelves(pat, lvls)
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

function fast_first_fit_decr!(
	sel :: AbstractArray{D},
	pat :: AbstractArray{D, 2},
	lvls :: Levels{D, S},
	p :: AbstractVector{P},
	l :: AbstractVector{S},
	w :: AbstractVector{S},
	t :: AbstractVector{D},
	n :: D, # number of distinct piece types (n >= maximum(t))
	L :: S, # length of the plate
	W :: S # width of the plate
) :: P where {D, S, P}
	c = lvls
	# assert explanation: all vectors should be the same length.
	#@assert isone(length(unique!(length.([p, l, w, t]))))
	# assert explanation: this is first-fit WIDTH-decreasing
	#@assert issorted(w; rev = true)
	#lw = Vector{S}() # levels width
	#ll = Vector{S}() # levels length
	rw = W # remaining width
	fill!(sel, zero(eltype(sel)))
	empty!(lvls)
	# The `pat` parameter is not reset because the code already
	# works on the assumption that it is a dirty buffer.

	bkv = zero(P) # best known value
	#sel = zeros(D, n) # if piece i is cut then sel[t[i]] is incremented
	#pat = Vector{Vector{D}}() # the pattern (always a 2-stage pattern)
	for i in eachindex(t) # for each piece
		assigned = false
		for j in 1:c.qt[] # for each already open level
			# assert reasoning: the levels are created with the width of previous
			# pieces, and the pieces are ordered by decreasing width, therefore
			# no piece may have more width than an already open level.
			#@assert lw[j] >= w[i]
			if c.ll[j] >= l[i] # checks if the level has enough length yet
				c.ll[j] -= l[i]
				bkv += p[i]
				sel[t[i]] += one(D)
				c.la[j] += 1 # increase the amount of pieces in that level
				pat[c.la[j], j] = t[i]
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
			c.qt[] += 1
			c.lw[c.qt[]] = w[i]
			c.ll[c.qt[]] = L - l[i]
			c.la[c.qt[]] = one(D)
			rw -= w[i]
			bkv += p[i]
			sel[t[i]] += one(D)
			pat[1, c.qt[]] = t[i]
		end
	end

	return bkv
end
end # module

