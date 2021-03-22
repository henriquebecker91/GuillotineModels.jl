"""
Collection of the methods used for generating the cuts and plates
that make up the PP-G2KP model.
"""
module Enumeration

# TODO: Consider removing flag ignore_2th_dim. Reasoning: nobody does that,
# even Furini only consider the pieces that fit a plate before discretizing it,
# even if ignoring the 2th dimension would save a lot of effort in the
# preprocessing nobody believes would be a valid trade-off to do so (the
# preprocessing phase participation on the total time is negligible), and it
# makes the code far more complex (because I have to have code that respects
# and ignores a basic assumption).

using DocStringExtensions
using UnPack: @unpack

using ...Utilities: SortedLinkedLW, switched_dims

export gen_cuts, ByproductPPG2KP
export should_extract_piece_from_plate

"""
    discretize(d::[D], l::[S], w::[S], L::S, W::S; kwargs...) :: [S]

!!! **Internal use.**

Given the pieces demand, length, and width (`d`, `l`, `w`), as well the plate
dimensions (`L`, `W`), return a vector of all linear combinations (constrained
by demand) of the length of pieces that fit the plate (the returned vector is
sorted by increasing value and no value is greater than `L`). To discretize the
width just swap `l` with `w` and `L` with `W`.

# Keyword arguments

All keyword arguments are of type `Bool` and have `false` as default.

* `only_single_pieces`: Instead of discretizing the length, return the
  piece lengths that exist only as the length of a piece and not as the
  combination of two or more smaller piece lengths.
* `ignore_W`: Do not use `W` to exclude pieces that do not fit the plate.
* `ignore_d`: Do not use `d`emand information to discretize (consider
  an unlimited amount of each piece type available).
"""
function discretize(
	d :: Vector{D}, l :: Vector{S}, w :: Vector{S}, L :: S, W :: S;
	only_single_pieces :: Bool = false, ignore_W :: Bool = false,
	ignore_d :: Bool = false
) where {D, S}
	# If two pieces have the same dimension they should be merged into a single
	# piece (summing their demands) for performance reasons.
	#@assert l == unique(l)
	# Each piece in 1:N has a demand and a dimension.
	@assert length(d) == length(l)
	N = length(d)
	# Remove items that cannot fit into L anyway.
	d_, l_, w_ = similar.((d, l, w))
	j = 0
	for i in 1:N
		if l[i] <= L && (ignore_W || w[i] <= W)
			j += 1
			d_[j] = d[i]
			l_[j] = l[i]
			w_[j] = w[i]
		end
	end
	resize!(d_, j)
	resize!(l_, j)
	resize!(w_, j)
	d, l, w = d_, l_, w_
	N = j
	iszero(N) && return S[]
	# marks: for each unit of the plate dimension, if there can be a cut there
	# (considering the pieces available) or not.
	marks = fill(false, L)
	# Mark the cuts of the first piece.
	!only_single_pieces && (marks[l[1]] = true)
	y = l[1] # y: used to iterate capacity, inherited from knapsack papers
	for _ = 2:(ignore_d ? L ÷ l[1] : min(d[1], L ÷ l[1]))
		marks[y += l[1]] = true
	end
	# If the lengths of the single pieces are not marked, the pairs are needed.
	if only_single_pieces
		for pii = 2:N, pij = 1:(pii - 1) # PIece `j` (as opposed to `i`)
			lij = l[pii] + l[pij]
			lij <= L && (marks[lij] = true)
		end
	end
	# Mark the cuts of all other pieces.
	for pii = 2:N # PIece Index
		li = l[pii] # length of i
		di = d[pii] # demand of i
		for y = (L - li):-1:1
			if marks[y]
				yrli = y + li # yrli: y + repeated lengths of i
				marks[yrli] = true
				for r = 2:(ignore_d ? L : di)
					yrli += li
					yrli > L && break
					marks[yrli] = true
				end
			end
		end
		# Mark cuts using just multiple copies of the same piece type (could be
		# done using a dummy cut at index zero with value true, but the array has
		# no index zero).
		y = li
		!only_single_pieces && (marks[y] = true)
		for _ = 2:(ignore_d ? L ÷ li : min(di, L ÷ li))
			marks[y += li] = true
		end
	end

	cuts = Vector{S}()
	sizehint!(cuts, L)
	if only_single_pieces
		# If we just want to check if the piece lengths cannot be obtained from
		# linear combinations of other pieces, we just iterate the unique piece
		# lengths in increasing order, and if they are not marked, then there is
		# no linear combination of other pieces that give the same length.
		usl = unique!(sort(l))
		for y in usl
			!marks[y] && push!(cuts, y)
		end
	else
		for (position, is_marked) in enumerate(marks)
			is_marked && push!(cuts, position)
		end
	end
	cuts
end

"""
    should_extract_piece_from_plate(pii, L, W, sllw, symm = 0x03) :: Bool

!!! **Internal use.**

Considering the enhanced PP-G2KP model, it is necessary to known if there
will be an 'extraction' variable representing the extraction of piece
with index `pii` from a plate with length `L` and width `L`.

A piece should be extracted from a plate if:
1. The piece fits inside the plate.
2. It is not possible to extract the piece in consideration together with any
   other piece (even other copy of the same piece) from the plate in
   consideration (for the cut orientations allowed by `symm` which by default
   are both vertical and horizontal).

# Arguments

1. `pii::D`: The index of the piece.
2. `L::S`: The length of the plate.
3. `W::S`: The width of the plate.
4. `sllw::SortedLinkedLW{D, S}`: The SortedLinkedLW struct for the pieces.
5. `symm::UInt8 = 0x03`: which cut orientations are allowed: `0x01` means
   'allow only horizontal cuts', `0x02` means 'allow only vertical cuts',
   and `0x03` means 'allow both kinds of cuts'.
"""
function should_extract_piece_from_plate(
	pii :: D, L :: S, W :: S, sllw :: SortedLinkedLW{D, S}, symm :: UInt8 = 0x03
) :: Bool where {D, S}
	li = sllw.l[pii]
	wi = sllw.w[pii]
	# If the piece does not fit the plate, then it is not even possible to
	# extract the piece from the plate.
	(li > L || wi > W) && return false

	n = length(sllw.l)

	# If we cut the plate by the length of the piece (a horizontal cut), does
	# this leaves space for another piece? If it does then we postpone to extract
	# the piece from a child of the current plate. Note that we can only do that
	# if the plate is marked to allow horizontal cuts or any cuts.
	# (The 'f' in the variables is from 'frontal'.)
	if symm == 1 || symm == 3
		fl = L - li
		fw = W

		for sli = 1:n # for the pieces in increase-or-stay order of their length
			# if the current piece does not fit the length none after it will fit
			sllw.sl[sli] > fl && break
			# if a piece fits the frontal child, then postpone extraction
			sllw.w[sllw.sli2pii[sli]] <= fw && return false
		end
	end

	if symm == 2 || symm == 3
		# The same as above but for the width of the piece (a vertical cut).
		# (The 'l' in the variables is from 'lateral'.)
		ll = L
		lw = W - wi

		# the same comments of the last loop are valid for this one
		for swi = 1:n
			sllw.sw[swi] > lw && break
			sllw.l[sllw.swi2pii[swi]] <= ll && return false
		end
	end

	# If no piece is small enough to make the plate-parameter be cut in a length
	# (or width) that leaves a second child large enough for extracting the
	# piece-parameter, then we cannot postpone and should allow extracting the
	# piece-parameter from the plate-parameter.
	return true
end

"""
    ub_num_pieces_fit(d, l, w, a, L, W, A, cutoff; ignore_2th_dim = false)

!!! **Internal use.**

Given a piece set (defined by `d`emand, `l`ength, `w`idth, and `a`rea),
and a plate (defined by `L`ength, `W`idth, and `A`rea), returns the `min`
between `cutoff` and an upper bound on the number of pieces that fit inside
the plate. If `ignore_2th_dim` is true, do not use `w` and `W` to filter
pieces that do not even fit the plate when alone.

In other words, `ub_num_pieces_fit(..., 6)` will return `6` if it is
possible (but not guaranteed) for `6` or more pieces to fit in the plate,
and will return `x` if `x < 6` is the greatest number of pieces that have
yet some possibility of being packed together; `x + 1` is already
guaranteed to be impossible.
"""
function ub_num_pieces_fit(
	d :: Vector{D},
	l :: Vector{S},
	w :: Vector{S},
	a :: Vector{P},
	L :: S,
	W :: S,
	A :: P,
	cutoff :: D;
	ignore_2th_dim = false
) :: D where {D, S, P}
	idxs = collect(1:length(a))
	sort!(idxs, by = i -> a[i])
	ub = zero(D)
	for i in idxs
		if l[i] <= L && (ignore_2th_dim || w[i] <= W)
			for _ = 1:d[i]
				A < a[i] && return ub
				A -= a[i]
				ub += 1
				ub > cutoff && return cutoff
			end
		end
	end
	return ub
end

"""
    no_chance_to_fit_6_piece(d, l, w, a, L, W, A; ignore_2th_dim = false)

Returns true if it is demonstrably impossible to fit six pieces into the plate;
and false if there is some possibility (but no guarantee) of fitting six pieces
in the plate.

If `ignore_2th_dim` is true, `w` and `W` are ignored by the algorithm.
"""
function no_chance_to_fit_6_piece(
	d :: Vector{D},
	l :: Vector{S},
	w :: Vector{S},
	a :: Vector{P},
	L :: S,
	W :: S,
	A :: P;
	ignore_2th_dim = false
) :: Bool where {D, S, P}
	ub = ub_num_pieces_fit(
		d, l, w, a, L, W, A, 6; ignore_2th_dim = ignore_2th_dim
	)
	return ub < 6
end

"""
    reduce2fit_usl(sl, sli2pii, w, L, W) :: [S]

!!! **Internal use.**

Given a plate (`L`, `W`), and the pieces sorted by length (`sl`, `sli2pii`,
`w`), return an increasing vector of unique lengths that pertain to a piece
that fits the plate. If two or more pieces share length and fit the plate, only
one copy of that length is included in the returned vector.

# Arguments

1. `sl::Vector{S}`: The piece lengths sorted by increase-or-stay order.
2. `sli2pii::Vector{D}`: If `sli2pii[i] == j` then `sl[i]` and `w[j]`
   correspond to the same item.
3. `w::Vector{S}`: The piece widths (in 'original' order).
4. `L::S`: The length of the plate.
5. `W::S`: The width of the plate.

"""
function reduce2fit_usl(
	sl :: Vector{S},
	sli2pii :: Vector{D},
	w :: Vector{S},
	L :: S,
	W :: S
) :: Vector{S} where {D, S, P}
	xs = empty(sl)
	for (i, x) in enumerate(sl)
		x > L && break # stop before the first lenght greater than L
		!isempty(xs) && last(xs) == x && continue # do not add duplicate lenghts
		w[sli2pii[i]] <= W && push!(xs, x) # add only if fits BOTH dimensions
	end

	xs
end

"""
    fits_at_least_one(sllw, L, W) :: Bool

!!! **Internal use.**

True if at least one piece in `sllw` fits inside a `L`x`W` plate;
false otherwise.
"""
function fits_at_least_one(
	sllw :: SortedLinkedLW{D, S}, L :: S, W :: S
) :: Bool where {D, S}
	for (i, v) in enumerate(sllw.sl)
		v > L && return false
		sllw.w[sllw.sli2pii[i]] <= W && return true
	end

	return false
end

"""
    filter_symm_pos(disc, dim, disc_copy = copy(disc)) -> disc_copy

!!! **Internal use.**

Implement the symmetry-breaking described in Section 2.1 of
10.1287/ijoc.2016.0710 (just after equation 10). Given a list of discretized
lengths (or widths) `disc` and the plate length (or width) `dim`, return a
copy of `disc` without all positions `i` in which ```disc[i] > div(dim/2)
\\land (\\exists. disc[j] == dim - disc[i])```.
"""
function filter_symm_pos(
	disc :: Vector{S},
	dim :: S,
	disc_copy = copy(disc) :: Vector{S}
) :: Vector{S} where {S}
	half_dim = dim ÷ 2
	for xy in disc
		(xy > half_dim || (xy == dim - xy)) && break
		symm_pos = dim - xy
		symm_idx = searchsortedlast(disc, symm_pos)
		!iszero(symm_idx) && disc[symm_idx] == symm_pos && (disc_copy[symm_idx] = 0)
	end
	return filter!(!iszero, disc_copy)
end

# From Furini2016 supplement: "if one (or more) of the flags of plate j has
# value -1, do not perform trim cuts on j in the flag orientation."
function filter_redundant_cuts!(
	# A list of cuts (parent plate, first child, second child). Only the
	# second child may be waste, and if it is, the index of the second child
	# is zero. All cuts share the same orientation (vertical or horizontal).
	nnn :: Vector{NTuple{3, P}},
	# Self Flag: corresponds to 'sh' or 'sv'. The orientation must match 'nnn'.
	sf :: Vector{Int8},
	# Father Flag: corresponds to 'fh' or 'fv'. Analogue to Self Flag.
	ff :: Vector{Int8}
) :: Vector{NTuple{3, P}} where {P}
	return filter!(nnn) do cut
		pp, _, sc = cut # parent plate, first child, second child
		# A cut should be kept if it is not a trim cut (!iszero(sc)) or,
		# in the case it is a trim cut, if none of its flags indicate a trim
		# cut was already made to obtain pp (check sf), or to obtain every possible
		# father plate of pp (check ff).
		!iszero(sc) || (sf[pp] > -1 && ff[pp] > -1)
	end
end

#=
Old documentation of the gen_cuts method. It was superseeded by the
ByproductPPG2KP documentation, but it is more complete than it and so
was kept here for now.

`pli_lwb`: List of all generated plates. The triple is: plate length, plate
	width, and plate bound (how many of the plate fit inside the original plate).
	The index/identifier of a plate is the same as its position in this list.
`hcuts`: If `(pp, cp1, cp2)` is present in this list, then the model has an
   horizontal cut that divides parent plate `pp` into child plates `cp1` and
   `cp2`. The three values are indexes in the list of all generated plates
   (except `cp2` may be zero, i.e., a trim/waste).
`vcuts`: The same as the second returned value, but for the vertical cuts.
`np`: A list of plate-piece pairs in which the plate may be sold as the
   respective piece. If `faithful2furini2016` is `true`, then this list
   has always the same size than the number of piece types and each
   plate that may be sold as a piece has the exact dimensions of the respective
   piece. If `faithful2furini2016` is `false` then this is the list of
   'extraction variables' and the plates may have dimensions slight
   different from their respective pieces, and each plate may be sold as
   one of many distinct pieces, and each piece may be extracted from many
   different plates.
=#
# TODO: Optional. Consider if changing the order of the variables can impact
# the solver performance. The only problem with this is that current
# ByproductPPG2KP forces all horizontal cuts to come before all vertical
# cuts. Changing the order of the constraints, however, needs the complete
# recomputation of ByproductPPG2KP.cuts and ByproductPPG2KP.np.
# Some criteria to be considered for the np variables (mostly):
# how much is the absolute profit obtained; how much relative area is wasted;
# how much is the profit of the a squared unit of the plate used; or group
# them by pii just for making it easier to the demand constraint.
"""
Collection of internal datastructures used to create the PP-G2KP model and
needed to assemble the solution given by the values of the model variables.

$TYPEDFIELDS

"""
struct ByproductPPG2KP{D, S, P}
	# pp (parent plate), fc (first child), sc (second child)
	"If `(pp, fc, sc)` is in `cuts`, then `pp` may be cut into `fc` and `sc`."
	cuts :: Vector{NTuple{3, P}}
	"If hybridize-with-restricted, has the extracted piece for each `cuts` elem."
	cut_extraction :: Vector{D}
	"Every cut before this position is horizontal, the rest are vertical."
	first_vertical_cut_idx :: P
	"If `(n, p)` is in `np` then plate `n` may be sold as piece `p`."
	np :: Vector{Tuple{P, D}} # TODO: rename np to pe (piece extraction)?
	"Indexed by a plate index, values are the plate length, width, and bound."
	pli_lwb :: Vector{Tuple{S, S, P}}
	"The demand of the pieces."
	d :: Vector{D}
	"The length of the pieces."
	l :: Vector{S}
	"The width of the pieces."
	w :: Vector{S}
	#p :: Vector{P} # for now, it is not necessary
	"The length of the original plate."
	L :: S
	"The width of the original plate."
	W :: S
end
Base.broadcastable(x :: ByproductPPG2KP{D, S, P}) where {D, S, P} = Ref(x)

# `_piece_only_positions` is similar to `guarantee_discretization!` but:
# 1. It is only necessary if hybridize-with-restricted is enabled.
# 2. It is eager instead of lazy (i.e., it computes the whole data
#    structure instead of filling the gaps as the enumeration goes).
#    The reasoning is that: (i) different from the discretization,
#    the number of distinct 'piece-only positions' is bounded by the
#    number of pieces (linear); (ii) it is almost certain that we will
#    need all of them.
# The return is a `Vector{Vector{S}}`, lets call it `pols`.
# `pows` outer Vector is indexed by the width available and gives the
# inner vector that is all the lengths that can only be obtained by
# single pieces (not by any combination) considering the length available.
# To obtain its counterpart, `pows`, indexed by the width available and
# returning the lengths only obtainable by single pieces, just switc
# `L` and `W` as well all fields of SortedLinkedLW for the call.
function _piece_only_positions(
	L :: S, W :: S, sllw :: SortedLinkedLW{D, S}, d; ignore_d = false
) :: Vector{Vector{S}} where {D, S}
	@unpack l, w = sllw
	pols = Vector{Vector{S}}(undef, W)
	empty_list = S[]
	# If less width is available than any piece the only valid value
	# is an empty list.
	for x in 1:(first(sllw.sw) - 1)
		pols[x] = empty_list
	end

	for (swi, available_w) in pairs(sllw.sw)
		not_is_last = swi < lastindex(sllw.sw)
		if not_is_last && sllw.sw[swi + 1] == available_w
			continue # i.e., only use the last in a block of repeated values
		end
		pol = discretize(
			d, l, w, L, available_w; only_single_pieces = true,
			ignore_d = ignore_d
		)
		before_next_w = not_is_last ? sllw.sw[swi + 1] - one(S) : W
		for x in available_w:before_next_w
			pols[x] = pol
		end
	end

	return pols
end

# Same as _piece_only_positions, but does the work for both dimensions
# instead of just the length dimension.
# The keyword arguments are passed to `_piece_only_positions` which
# pass them to `discretize`.
function piece_only_positions(
	L :: S, W :: S, sllw :: SortedLinkedLW{D, S}, d; kwargs...
) :: NTuple{2, Vector{Vector{S}}} where {D, S}
	pols = _piece_only_positions(L, W, sllw, d; kwargs...)
	rev_sllw = switched_dims(sllw)
	pows = _piece_only_positions(W, L, rev_sllw, d; kwargs...)

	return pols, pows
end

function guarantee_discretization!(
	dls, dws, pll :: S, plw :: S, L :: S, W :: S, d, l, w,
	round2disc, ignore_2th_dim, ignore_d, level :: UInt8 = UInt8(1)
) :: Tuple{S, S} where {S}
	# If our understanding of the discretization is correct, it is impossible
	# for this method to recurse more than one time in a row.
	@assert level == 1 || level == 2
	(iszero(pll) || iszero(plw)) && return zero(S), zero(S)

	# First of all, let us get the two discretizations, memoized or not.
	dls_plw = if isassigned(dls, plw)
		dls[plw]
	else
		discretize(d, l, w, L, plw; ignore_W = ignore_2th_dim, ignore_d = ignore_d)
	end
	dws_pll = if isassigned(dws, pll)
		dws[pll]
	else
		discretize(d, w, l, W, pll; ignore_W = ignore_2th_dim, ignore_d = ignore_d)
	end

	# The code for fully initializing `dls[plw]` needs partially
	# initializing `dws[L]`, if we also wanted to fully initialize
	# `dws[L]` then we need to save the information that it is not
	# fully initialized (to remember initializing it after).
	dws_L_full_init_pending = pll == L && !isassigned(dws, L)

	# The following two `if`s do essentially the same thing but one is for `dlw`
	# and the other for `dlw`. They find a range of the second (non-discretized)
	# dimension for which the discretization (of the first dimension) is valid.
	# Then assign the same discretization to all positions in the range, possibly
	# saving the effort for new calls (an eager memoization).
	if !isassigned(dls, plw)
		dws_L_full_init_pending && (dws[L] = dws_pll)
		idx_lte_plw_in_dws_L = searchsortedlast(dws[L], plw)
		if iszero(idx_lte_plw_in_dws_L)
			# `searchsortedlast` returns zero if the first value in `dws[L]` is
			# already greater than `plw`, i.e., there is no piece or piece set
			# that fit in so little width (the minimum width, or first
			# discretized position, is after `plw`).
			@assert isempty(dls_plw) # can only be an empty discretization set
			# Assign the empty set to these positions so we to not lose time
			# recomputing them. Those are the first spaces no piece fits.
			!isempty(dws[L]) && (dls[1:(first(dws[L]) - 1)] .= (dls_plw,))
		elseif idx_lte_plw_in_dws_L == length(dws[L])
			# If `idx_lte_plw_in_dws_L` is the last discretization in `dws[L]`
			# this means that from there to `W` all positions can only pack the
			# same pieces/'piece sets'.
			first_valid_width = dws[L][idx_lte_plw_in_dws_L]
			dls[first_valid_width:W] .= (dls_plw,)
		else
			# If `idx_lte_plw_in_dws_L` is nor zero, nor the last discretization,
			# then `dls_plw` is valid for the range `dws[L][idx_lte_plw_in_dws_L]`
			# to `dws[L][idx_lte_plw_in_dws_L + 1] - 1` (i.e., one width unit
			# before the next discretized width position).
			first_valid_width = dws[L][idx_lte_plw_in_dws_L]
			last_valid_width = dws[L][idx_lte_plw_in_dws_L + 1] - 1
			dls[first_valid_width:last_valid_width] .= (dls_plw,)
		end
	end
	if !isassigned(dws, pll) || dws_L_full_init_pending
		@assert plw != W || isassigned(dls, W)
		idx_lte_pll_in_dls_W = searchsortedlast(dls[W], pll)
		# See comments of the similar `if` block above.
		if iszero(idx_lte_pll_in_dls_W)
			@assert isempty(dws_pll)
			!isempty(dls[W]) && (dws[1:(first(dls[W]) - 1)] .= (dws_pll,))
		elseif idx_lte_pll_in_dls_W == length(dls[W])
			first_valid_length = dls[W][idx_lte_pll_in_dls_W]
			dws[first_valid_length:L] .= (dws_pll,)
		else
			first_valid_length = dls[W][idx_lte_pll_in_dls_W]
			last_valid_length = dls[W][idx_lte_pll_in_dls_W + 1] - 1
			dws[first_valid_length:last_valid_length] .= (dws_pll,)
		end
	end

	if round2disc
		r2d_pll_idx = searchsortedlast(dls[plw], pll)
		r2d_plw_idx = searchsortedlast(dws[pll], plw)
		r2d_pll = iszero(r2d_pll_idx) ? zero(S) : dls[plw][r2d_pll_idx]
		r2d_plw = iszero(r2d_plw_idx) ? zero(S) : dws[pll][r2d_plw_idx]
		if iszero(r2d_pll) || iszero(r2d_plw)
			r2d_pll = r2d_plw = zero(S)
		end
		#println("Before recursion")
		#@show plw, pll, r2d_plw, r2d_pll
		r2d_not_guaranteed = (!iszero(r2d_pll) && !isassigned(dws, r2d_pll)) ||
			(!iszero(r2d_plw) && !isassigned(dls, r2d_plw))
		if r2d_not_guaranteed
			r2d_pll, r2d_plw = guarantee_discretization!(
				dls, dws, r2d_pll, r2d_plw, L, W, d, l, w,
				round2disc, ignore_2th_dim, ignore_d, level + one(level)
			)
		end
		@assert iszero(r2d_pll) || isassigned(dws, r2d_pll)
		@assert iszero(r2d_plw) || isassigned(dls, r2d_plw)
		return r2d_pll, r2d_plw
	else
		return pll, plw
	end
end

# mirror plate
_mp(::Val{true}, l, w) = l > w ? (w, l) : (l, w)
_mp(::Val{false}, l, w) = l, w

macro mp(v, l, w)
	v_, l_, w_ = esc(v), esc(l), esc(w)#gensym(v), gensym(l), gensym(w)# v, l, w #
	return :(_mp($v_, $l_, $w_)...)
end

function do_cut!(
	is_vertical :: Bool, # If the cut is vertical (i.e., on a width position)
	n :: P, # The last plate id, incremented with each plate created and returned.
	m :: Val{M}, # If the plates are mirrored (i.e., length always smaller).
	fcl :: S, # The length of the first child (if round2disc then it may shrink).
	fcw :: S, # The width of the first child (if round2disc then it may shrink).
	scl :: S, # The length of the second child (if round2disc then it may shrink).
	scw :: S, # The width of the second child (if round2disc then it may shrink).
	next_idx :: P, # The index of next and sfhv for the current parent plate.
	dls :: Vector{Vector{S}}, # Discretized lengths set for each width.
	dws :: Vector{Vector{S}}, # Discretized widths set for each length.
	pols :: Vector{Vector{S}}, # Piece-Only Length Sets (by available width).
	pows :: Vector{Vector{S}}, # Piece-Only Widths Sets (by available length).
	L :: S, # The original plate length.
	W :: S, # The original plate width.
	d :: Vector{D}, # The pieces demand (needed for guaranteeing discretization).
	l :: Vector{S}, # The pieces length.
	w :: Vector{S}, # The pieces width.
	sllw :: SortedLinkedLW{D, S}, # Used for `fits_at_least_one`
	sfhv :: NTuple{4, Vector{Int8}}, # Used for the Redundant-Cut optimization.
	next :: Vector{Tuple{S, S, P}}, # New plates are to be pushed here.
	plis :: Array{P, 2}, # plis[pll, plw] == pli (if plate already known)
	nnn :: Vector{NTuple{3, P}}, # Either vnnn or hnnn.
	dce :: Vector{D}, # Either vdce or hdce.
	round2disc :: Bool, # Passed to guarantee_discretization!
	ignore_2th_dim :: Bool, # Passed to guarantee_discretization!
	ignore_d :: Bool, # Passed to guarantee_discretization!
	faithful2furini2016 :: Bool,
	hybridize_with_restricted :: Bool,
	aggressive_hybridization :: Bool,
	no_redundant_cut :: Bool,
) :: Tuple{P, P, P, P} where {D, S, P, M}
	aggressive_hybridization_extra_vars = zero(P)
	# fcl: first child length
	# fcw: first child width
	fcl, fcw = guarantee_discretization!(
		dls, dws, fcl, fcw, L, W, d, l, w,
		round2disc, ignore_2th_dim, ignore_d
	)
	# scl: second child length
	# scw: second child width
	scl, scw = guarantee_discretization!(
		dls, dws, scl, scw, L, W, d, l, w,
		round2disc, ignore_2th_dim, ignore_d
	)

	# dcpiis: Double-Cut PIece IndexeS. Has a single zero if no hybridization
	#   has happened, has one or more piece ids if hybridizations have happened.
	# fcdpars: if multiple hybridizations happen, then multiple `fcdpar` exist,
	#   so after the `if hybridize_with_restricted` block below, we loop over
	#   `fcdpars`. If there was no hybridization, then the only `fcdpar` is
	#   pushed into the vector; if there is one or more hybridizations, they
	#   are pushed into the `fcdpar`.
	dcpiis = D[]
	fcdpars = S[]
	# pop -- Piece-Only Position (always perpendicular to cut)
	# fcdper -- First Child Dimension PERpendicular to cut
	# fcdpar -- First Child Dimension PARallel to cut
	# pii2dper -- PIece Index TO Dimension PERpendicular to cut
	# pii2dpar -- PIece Index TO Dimension PARallel to cut
	# sddper -- Sorted Discretized Dimensions PERpendicular to cut
	# sddpar -- Sorted Discretized Dimensions PARallel to cut
	# sdd2pii -- Sorted Discretized Dimensions TO PIece Index (perpendicular)
	# First rename the dimensions to be cut-oriented.
	fcdper, fcdpar = is_vertical ? (fcw, fcl) : (fcl, fcw)
	if hybridize_with_restricted
		# Now rename all the other relevant data structures to be cut-oriented.
		pop, pii2dper, pii2dpar, sddper, sddpar, sdd2pii = if is_vertical
			(ignore_2th_dim ? pows[L] : pows[fcl]), sllw.w, sllw.l,
				sllw.sw, sllw.sl, sllw.swi2pii
		else
			(ignore_2th_dim ? pols[W] : pols[fcw]), sllw.l, sllw.w,
				sllw.sl, sllw.sw, sllw.sli2pii
		end
		# If `fcdper` is in `pop`, then we have one or more hybridizations,
		# otherwise we have no hybridization.
		fcdper_pop_idx = searchsorted(pop, fcdper)
		if isempty(fcdper_pop_idx) # no hybridization
			push!(dcpiis, 0)
			push!(fcdpars, fcdpar)
		else # one or more hybridizations
			spiis = searchsorted(sddper, fcdper)
			# If there is a single hybridization possibility, then we always
			# hybridize. However, if there are multiple possibilities, we only
			# hybridize if aggressive_hybridization is enabled, to avoid
			# increasing the model size too much.
			if isone(length(spiis)) || aggressive_hybridization
				for spii in spiis
					pii = sdd2pii[spii]
					if pii2dpar[pii] <= fcdpar
						push!(dcpiis, convert(D, pii))
						push!(fcdpars, fcdpar - pii2dpar[pii])
					end
				end
				aggressive_hybridization_extra_vars += length(spiis) - 1
			else # multiple hybridization candidates but !aggressive_hybridization
				# Abstain from any hybridization to keep optimality and model size.
				push!(dcpiis, 0)
				push!(fcdpars, fcdpar)
			end
			@assert !isempty(dcpiis)
			@assert !isempty(fcdpars)
		end
	else
		# If !hybridize_with_restricted, we fill dcpiis and fcdpars with
		# sentinel values, that are the correct if no hybridization has ocurred.
		push!(dcpiis, 0)
		push!(fcdpars, fcdpar)
	end # if hybridize_with_restricted

	qt_first_born_kills :: P = zero(P)
	for (dcpii, fcdpar) in zip(dcpiis, fcdpars)
		fcl, fcw = is_vertical ? (fcdpar, fcdper) : (fcdper, fcdpar)

		# If `dcpii` is not zero, then this is an hybridization, and the plate
		# may be further reduced by round2disc.
		if !iszero(dcpii) && round2disc
			fcl, fcw = guarantee_discretization!(
				dls, dws, fcl, fcw, L, W, d, l, w,
				round2disc, ignore_2th_dim, ignore_d
			)
		end

		n, killed_first_born = do_cut!(is_vertical,
			n, m, fcl, fcw, scl, scw, next_idx, dls, dws,
			sllw, sfhv, next, plis, nnn, faithful2furini2016,
			hybridize_with_restricted, no_redundant_cut,
		)
		qt_first_born_kills += killed_first_born

		# Finally, if hybridize_with_restricted is enabled, and the current cut
		# was transformed by this optimization in a double cut, we need to save
		# the extracted piece in dce (otherwise zero is our sentinel).
		hybridize_with_restricted && push!(dce, dcpii)
	end

	qt_hybridizations :: P = zero(P)
	if !iszero(first(dcpiis))
		qt_hybridizations = convert(P, length(dcpiis))
	end

	return (n, qt_hybridizations, qt_first_born_kills,
		aggressive_hybridization_extra_vars)
end

function do_cut!(
	is_vertical :: Bool, # If the cut is vertical (i.e., on a width position)
	n :: P, # The last plate id, incremented with each plate created and returned.
	m :: Val{M}, # If the plates are mirrored (i.e., length always smaller).
	fcl :: S, # The length of the first child (if round2disc then it may shrink).
	fcw :: S, # The width of the first child (if round2disc then it may shrink).
	scl :: S, # The length of the second child (if round2disc then it may shrink).
	scw :: S, # The width of the second child (if round2disc then it may shrink).
	next_idx :: P, # The index of next and sfhv for the current parent plate.
	dls :: Vector{Vector{S}}, # Discretized lengths set for each width.
	dws :: Vector{Vector{S}}, # Discretized widths set for each length.
	sllw :: SortedLinkedLW{D, S}, # Used for `fits_at_least_one`
	sfhv :: NTuple{4, Vector{Int8}}, # Used for the Redundant-Cut optimization.
	next :: Vector{Tuple{S, S, P}}, # New plates are to be pushed here.
	plis :: Array{P, 2}, # plis[pll, plw] == pli (if plate already known)
	nnn :: Vector{NTuple{3, P}}, # Either vnnn or hnnn.
	faithful2furini2016 :: Bool,
	hybridize_with_restricted :: Bool,
	no_redundant_cut :: Bool,
) :: Tuple{P, Bool} where {D, S, P, M}
	# Put some info readily available.
	pll, plw, pli = next[next_idx] # PLate Length, Width, and Index
	if !no_redundant_cut
		sh, sv, fh, fv = sfhv
		sh_j, sv_j, fh_j, fv_j = getindex.(sfhv, next_idx)
	end
	# If the preprocessing is faithful2furini2016, the plates are cut until
	# they have the same size as pieces and, consequently, there exist the
	# concept of trim cut (i.e., a cut in which the second child is waste).
	trim_cut = faithful2furini2016 && !fits_at_least_one(sllw, scl, scw)
	# The trim_cut flag is used below outside of 'if faithful2furini2016'
	# because a true value implicates that 'faithful2furini2016 == true'.

	# For now, the only reason for a first child to be unable to pack
	# any piece is the hybridize_with_restricted optimization.
	killed_first_born = false
	if hybridize_with_restricted && !fits_at_least_one(sllw, fcl, fcw)
		killed_first_born = true
	elseif iszero(plis[@mp m fcl fcw]) # If the first child does not yet exist.
		push!(next, (fcl, fcw, n += 1)) # Create the first child.
		if !no_redundant_cut
			# For now, we do not know if hybridize_with_restricted is
			# entirely compatible with Redundant-Cut, so we only consider
			# a trim_cut if it does not clash with hybridize_with_restricted.
			if trim_cut && !hybridize_with_restricted
				# From Furini2016 supplement: "if a NEW plate j_1 ∈ J is obtained
				# from j through a trim cut with orientation v:"
				if is_vertical
					fchild_sfhv = (0, -1, sh_j, sv_j)
				else
					fchild_sfhv = (-1, 0, sh_j, sv_j)
				end
			else
				# If a plate (NEW or existing) j_1 is obtained from j without a trim
				# cut: set all flags of j 1 to 1.
				fchild_sfhv = (1, 1, 1, 1)
			end
			push!.(sfhv, fchild_sfhv)
		end
		plis[@mp m fcl fcw] = n # Mark plate existence.
	elseif !no_redundant_cut
		# If plate already exists, and we are implementing Redundant-Cut.
		# Here we do not care for hybridize_with_restricted because
		# setting positions to positive (i.e., `1`) is disabling
		# Redundant-Cut optimizations not enabling them.
		if trim_cut
			# From Furini2016 supplement: "if an existing plate j_1 ∈ J is
			# obtained from j through a trim cut:"
			sh_j > -1 && (fh[plis[@mp m fcl fcw]] = 1)
			sv_j > -1 && (fv[plis[@mp m fcl fcw]] = 1)
		else
			# From Furini2016 supplement: "if a plate (new or EXISTING) j_1 is
			# obtained from j without a trim cut: set all flags of j_1 to 1."
			setindex!.(sfhv, (1, 1, 1, 1), plis[@mp m fcl fcw])
		end
	end

	# If the second child is not waste (!trim_cut), nor was already
	# generated, then it is a new plate obtained without a trim cut.
	if !trim_cut && iszero(plis[@mp m scl scw])
		# Assert meaning: if the second child is not waste, then it must
		# have an associated size.
		@assert !iszero(scw)
		# Save the plate to be processed later, and mark its existence.
		push!(next, (scl, scw, n += 1))
		plis[@mp m scl scw] = n
		# From Furini2016 supplement: "if a plate (NEW or existing) j 1 is
		# obtained from j without a trim cut: set all flags of j_1 to 1."
		!no_redundant_cut && push!.(sfhv, (1, 1, 1, 1))
	end
	# Add the cut to the cut list. Check if the second child plate is waste
	# before trying to get its plate index.
	push!(nnn, (pli,
		killed_first_born ? 0 : plis[@mp m fcl fcw],
		trim_cut ? 0 : plis[@mp m scl scw]
	))

	return n, killed_first_born
end

"""
    gen_cuts(::Type{P}, d, sllw, L, W; [kwargs])

!!! **Internal use.**

The main method responsible for generating all the cuts used to create both the
original PP-G2KP model as its enhanced version (in this case, it generates the
piece extractions too).

# Arguments

## Positional and required arguments

1. `::Type{P}`: The type of the piece profits (and area).
2. `d::Vector{D}`: The pieces demand.
3. `sllw::SortedLinkedLW{D, S}`: The SortedLinkedLW containing the length and
   width of the pieces.
4. `L::S`: The plate length.
5. `W::S`: The plate width.

## Keyword arguments

All keyword arguments are of type `Bool`.

* `ignore_2th_dim`: Default: false. Ignore the dimension not being discretized
  during discretization.
* `ignore_d`: Default: false. Ignore the demand information during
  discretization.
* `round2disc`: Default: true. Round the size of the second child of a
  cut to a discretized position.
* `faithful2furini2016`: Default: false. Tries to be the most faithful possible
  to the description in 10.1287/ijoc.2016.0710.
* `no_redundant_cut`: Default: false. Disables the Redundant-Cut reduction
  described in 10.1287/ijoc.2016.0710.
* `no_cut_position`: Default: false. Disables the Cut-Position reduction
  described in 10.1287/ijoc.2016.0710.
* `no_furini_symmbreak`: Default: false. Ignored if `faithful2furini2016`
  is `false`. Disables the symmetry-breaking used in 10.1287/ijoc.2016.0710
  (consequently all discretized positions are used, none is removed).
"""
function gen_cuts(
	::Type{P}, d :: Vector{D}, sllw :: SortedLinkedLW{D, S}, L :: S, W :: S,
	m :: Val{M} = Val(false);
	ignore_2th_dim = false, ignore_d = false, round2disc = true,
	hybridize_with_restricted = false,
	aggressive_hybridization = false,
	faithful2furini2016 = false,
	no_redundant_cut = false, no_cut_position = false,
	no_furini_symmbreak = false, quiet = false, verbose = false
) :: ByproductPPG2KP{D, S, P} where {D, S, P, M}
	mirror_plates = M # the Bool inside the Val
	# If both verbose and quiet are passed, then quiet wins.
	verbose = verbose & !quiet
	!quiet && faithful2furini2016 && round2disc && @warn(
		"Enabling both faithful2furini2016 and round2disc is allowed, but you" *
		" are not being entirely faithful to Furini2016 if you do so."
	)
	!faithful2furini2016 && no_redundant_cut && !quiet && @warn(
		"The Redundant-Cut is only used when faithful2furini2016" *
		" is enabled. As flag faithful2furini2016 is not enabled, flag" *
		" no-redundant-cut has no effect."
	)
	if !no_redundant_cut && mirror_plates
		!quiet && @warn(
			"Flag mirror-plates has disabled the Redundant-Cut reduction."
		)
		no_redundant_cut = true
	end

	l = sllw.l
	w = sllw.w
	@assert length(d) == length(l)
	@assert length(d) == length(w)
	# unique versions
	usl = unique(sllw.sl)
	usw = unique(sllw.sw)
	a = l .* w
	A = L * W
	max_piece_type = convert(D, length(l))
	#max_num_plates = convert(P, L) * convert(P, W)
	min_pil = minimum(l)
	min_piw = minimum(w)

	# If (n1, n2, n3) is in nnn, then n1 may be partitioned in n2 and n3.
	# Not every possible partition is present. If the model is trying to be
	# the most faithful to Furini, then redundant-cuts is implemented. If the
	# model is not restricted to Furini's description, then it may cut even more.
	hnnn = Vector{NTuple{3, P}}()
	vnnn = Vector{NTuple{3, P}}()
	# dce (double cut extractions): if !hybridize_with_restricted these are kept
	# as empty Vectors; otherwise, some cuts are double cuts, that generate
	# up to two child plates as usual, but also a piece extraction. If vdce[i]
	# is zero, then hnnn[i] is not a double cut, otherwise hnnn[i] is a double
	# cut and hdci[i] indicate the piece extracted. Analogue for vnnn and vdce.
	hdce = Vector{D}()
	vdce = Vector{D}()
	# The np vector has the pairs plate-piece for which the respective plate is
	# allowed to be sold as the respective piece.
	# If the model is faithful to Furini, it just has a single plate for each
	# piece type and the plate has the exact size of the piece.
	# Otherwise, if the piece fits the plate size and the plate does not allow
	# any other piece by its side (i.e., there is no second piece that can be
	# placed in the same plate vertically or horizontally) then (plate, piece) is
	# in np.
	np = Vector{Tuple{P, D}}()
	# The list of plates attributes: plate length, plate width, and plate bound.
	# The plate index is the same as the index in pli_lwb.
	pli_lwb = Vector{Tuple{S, S, P}}()
	# plis: matrix of the plate dimensions in which zero means "never seen that
	# plate before" and nonzero means "this nonzero number is the plate index".
	if mirror_plates
		plis = zeros(P, max(L, W), max(L, W))
	else
		plis = zeros(P, L, W)
	end
	plis[@mp m L W] = one(P)
	# next: plates already indexed but not yet processed, starts with (L, W, 1).
	# Storing the index as the third value is not necessary (as it could be
	# queried from plis) but this is probably more efficient this way.
	next = Vector{Tuple{S, S, P}}()
	next_idx = one(P)
	#sizehint!(next, max_num_plates)
	push!(next, (L, W, one(P)))
	sfhv = (Vector{Int8}(), Vector{Int8}(), Vector{Int8}(), Vector{Int8}())
	if !no_redundant_cut
		# If the preprocessing is faithful2furini2016 and Redundant-Cut is
		# enabled, then four auxiliary trim cut flag vectors are necessary.
		sh, sv, fh, fv = sfhv
		push!.(sfhv, (1, 1, 1, 1))
	end
	# n: The amount of plates (the index of the highest plate type).
	n = one(P) # there is already the original plate
	# Memoized discretizations. The discretized lengths for every plate width.
	# The discretized widths for every plate length.
	dls = Vector{Vector{S}}(undef, W)
	dws = Vector{Vector{S}}(undef, L)
	# piece only length/width -- for use with hybridize_with_restricted,
	# has the length and widths that only appear on pieces and not on
	# piece combinations. We start simple, without the optimization below.
	# TODO: `pol` and `pow` could be a set of discretizations like `dls`/`dws`
	# because in smaller plates some piece length/width may stop being
	# attainable by combinations of smaller plates (e.g., a set of pieces
	# matched another piece length, but in the smaller plate some piece of
	# the set does not fit because of the width, while the piece that had
	# the length matched yet fits at that width).
	if hybridize_with_restricted
		pols, pows = piece_only_positions(L, W, sllw, d; ignore_d)
	else
		# If `hybridize_with_restricted` these will never be used, but
		# it is practical if they exist and are the right type to be passed
		# as arguments (and type-stable).
		pols = pows = Vector{S}[]
	end

	# There are two discretizations for each dimension: dl1/dw1 and dl2/dw2.
	# They are only different if the Cut-Position reduction is being used. In
	# this case, dl1/dw1, which is used for cutting the plate, may be a
	# restricted cut set (i.e., it does not have piece size combinations, just
	# single piece sizes). The dl2/dw2 is always the complete discretization,
	# and may be needed (even if Cut-Position is used) for reducing the size of
	# the second child to the last discretized position (if round2disc is
	# enabled).
	guarantee_discretization!(
		dls, dws, L, W, L, W, d, l, w, round2disc, ignore_2th_dim, ignore_d
	)
	hybridizations = 0
	first_borns_killed_by_hybridization = 0
	aggressive_hybridization_extra_vars = 0
	while next_idx <= length(next)
		#@show next_idx
		#@show length(next)
		@assert (no_redundant_cut || [length(next)] == unique(length.(sfhv)))

		pll, plw, pli = next[next_idx] # PLate Length, Width, and Index
		if !faithful2furini2016
			for pii in 1:max_piece_type # pii: PIece Index
				if should_extract_piece_from_plate(pii, pll, plw, sllw)
					push!(np, (pli, pii))
				end
			end
		end

		plb = (L ÷ pll) * (W ÷ plw) # PLate Bound
		# It is not necessary to store the plate id in pli_lwb because they are added
		# in order (a queue is used), so the array index is the plate index.
		push!(pli_lwb, (pll, plw, plb))
		if ignore_2th_dim
			dl1 = dl2 = dls[W]
			dw1 = dw2 = dws[L]
		else
			#@show plw
			@assert isassigned(dls, plw)
			#@show pll
			@assert isassigned(dws, pll)
			#=
			pll, plw = guarantee_discretization!(
				dls, dws, pll, plw, L, W, d, l, w, round2disc, ignore_2th_dim, ignore_d
			)
			=#
			dl1 = dl2 = dls[plw]
			dw1 = dw2 = dws[pll]
		end

		if !no_cut_position
			if no_chance_to_fit_6_piece(
				d, l, w, a, pll, plw, pll * plw; ignore_2th_dim = ignore_2th_dim
			)
				if ignore_2th_dim
					dl1 = usl
					dw1 = usw
				else
					pll_ = faithful2furini2016 ? pll : pll ÷ 2
					plw_ = faithful2furini2016 ? plw : plw ÷ 2
					dl1 = reduce2fit_usl(sllw.sl, sllw.sli2pii, w, pll_, plw)
					@assert dl1 ⊆ dl2
					dw1 = reduce2fit_usl(sllw.sw, sllw.swi2pii, l, plw_, pll)
					@assert dw1 ⊆ dw2
				end
				#= Debug. Showing something is not an error.
				if dl1 != dl2[1:length(dl1)] || dw1 != dw2[1:length(dw1)]
					fpr_count += 1
					@show pll
					@show plw
					if dl1 != dl2[1:length(dl1)]
						@show dl1
						@show dl2
					end
					if dw1 != dw2[1:length(dw1)]
						@show dw1
						@show dw2
					end
				end
				=#
			end
		end

		# If faithful2furini2016 flag is enabled, then the whole discretization is
		# iterated, except by the cuts on the second half of the plate that have an
		# exact symmetry on the fist half (same distance from plate midpoint). The
		# calls below create a new vector without such symmetric points.
		if faithful2furini2016 && !no_furini_symmbreak
			dl1 = filter_symm_pos(dl1, pll)
			dw1 = filter_symm_pos(dw1, plw)
		end

		for x in dw1
			# See comment above filter_symm_pos calls.
			if faithful2furini2016
				x >= plw && break
			else
				x > plw ÷ 2 && break
			end

			@assert pll >= min_pil
			@assert x >= min_piw
			@assert faithful2furini2016 || plw - x >= min_piw

			fcl, fcw, scl, scw = pll, x, pll, plw - x
			n, qt_hybridizations, qt_first_born_kills, agg_hyb_extra = do_cut!(true,
				n, m, fcl, fcw, scl, scw, next_idx, dls, dws, pols, pows,
				L, W, d, l, w, sllw, sfhv, next, plis,
				vnnn, # Either vnnn or hnnn, i.e., the correct orientation of the two
				vdce, # Either vdce or hdce, i.e., the correct orientation of the two
				round2disc, ignore_2th_dim, ignore_d,
				faithful2furini2016, hybridize_with_restricted,
				aggressive_hybridization, no_redundant_cut
			)
			hybridizations += qt_hybridizations
			first_borns_killed_by_hybridization += qt_first_born_kills
			aggressive_hybridization_extra_vars += agg_hyb_extra
		end

		for y in dl1
			# See comment above filter_symm_pos calls.
			if faithful2furini2016
				y >= pll && break
			else
				y > pll ÷ 2 && break
			end

			@assert plw >= min_piw
			@assert y >= min_pil
			@assert faithful2furini2016 || pll - y >= min_pil

			fcl, fcw, scl, scw = y, plw, pll -y, plw
			n, qt_hybridizations, qt_first_born_kills, agg_hyb_extra = do_cut!(false,
				n, m, fcl, fcw, scl, scw, next_idx, dls, dws, pols, pows,
				L, W, d, l, w, sllw, sfhv, next, plis,
				hnnn, # Either vnnn or hnnn, i.e., the correct orientation of the two
				hdce, # Either vdce or hdce, i.e., the correct orientation of the two
				round2disc, ignore_2th_dim, ignore_d,
				faithful2furini2016, hybridize_with_restricted,
				aggressive_hybridization, no_redundant_cut
			)

			hybridizations += qt_hybridizations
			first_borns_killed_by_hybridization += qt_first_born_kills
			aggressive_hybridization_extra_vars += agg_hyb_extra
		end

		# Goes to the next unprocessed plate.
		next_idx += 1
	end

	if faithful2furini2016
		# Apply Redundant-Cut.
		if !no_redundant_cut
			# red_cut_flags = collect(zip(sh, sv, fh, fv))
			# println("red_cut_flags begin")
			# show(stdout, "text/plain", red_cut_flags)
			# println()
			# println("red_cut_flags end")
			qt_cuts_before_redundant_cut = length(hnnn) + length(vnnn)
			filter_redundant_cuts!(hnnn, sh, fh)
			filter_redundant_cuts!(vnnn, sv, fv)
			qt_cuts_after_redundant_cut = length(hnnn) + length(vnnn)
			qt_cuts_removed_by_redundant_cut = qt_cuts_before_redundant_cut -
				qt_cuts_after_redundant_cut
			if verbose
				print(
					"""
					qt_cuts_before_redundant_cut = $qt_cuts_before_redundant_cut
					qt_cuts_after_redundant_cut = $qt_cuts_after_redundant_cut
					qt_cuts_removed_by_redundant_cut = $qt_cuts_removed_by_redundant_cut
					"""
				)
			end
		end

		# Create conversion table of the index of the plate with same size as
		# a piece to the corresponding piece index.
		for pii in 1:max_piece_type # pii: PIece Index
			# Assert meaning: for every distinct piece dimensions there should be
			# a plate with the exact same dimensions.
			#@assert !iszero(get(plis, (l[pii], w[pii]), 0))
			#push!(np, (plis[(l[pii], w[pii])], pii))
			@assert !iszero(plis[@mp m l[pii] w[pii]])
			push!(np, (plis[@mp m l[pii] w[pii]], pii))
		end
	end

	# The two commented 'for' loops below print all variable indexes. They were
	# used to check the variable generated against the Martin's code. They were
	# kept commented given the necessity arises again, but may be removed in the
	# future (or a flag may be created to output them).
	#=
	for hcut in hnnn
		pp, fc, _ = hcut
		ppl, ppw, _ = pli_lwb[pp]
		q = pli_lwb[fc][1]
		println("$(ppl) $(ppw) 1 $(q)")
	end
	for vcut in vnnn
		pp, fc, _ = vcut
		ppl, ppw, _ = pli_lwb[pp]
		q = pli_lwb[fc][2]
		println("$(ppl) $(ppw) 0 $(q)")
	end
	=#

	if verbose
		@show hybridizations
		@show first_borns_killed_by_hybridization
		@show aggressive_hybridization_extra_vars
	end

	first_vertical_cut_idx = length(hnnn) + 1
	return ByproductPPG2KP(
		append!(hnnn, vnnn), append!(hdce, vdce),
		first_vertical_cut_idx,
		np, pli_lwb, deepcopy(d), deepcopy(l), deepcopy(w),
		deepcopy(L), deepcopy(W)
	)
end

end # module

