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

import ...Utilities.SortedLinkedLW

export gen_cuts, ByproductPPG2KP
export should_extract_piece_from_plate

"""
    discretize(d::[D], l::[S], w::[S], L::S, W::S; kwords...) :: [S]

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

"""
    gen_cuts(::Type{P}, d, sllw, L, W; [kword_args])

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
	::Type{P}, d :: Vector{D}, sllw :: SortedLinkedLW{D, S}, L :: S, W :: S;
	ignore_2th_dim = false, ignore_d = false, round2disc = true,
	faithful2furini2016 = false,
	no_redundant_cut = false, no_cut_position = false,
	no_furini_symmbreak = false, quiet = false, verbose = false
) :: ByproductPPG2KP{D, S, P} where {D, S, P}
	# If both verbose and quiet are passed, then quiet wins.
	verbose = verbose & !quiet
	!quiet && faithful2furini2016 && round2disc && @warn(
		"Enabling both faithful2furini2016 and round2disc is allowed, but you" *
		" are not being entirely faithful to Furini2016 if you do so."
	)
	if !faithful2furini2016 && no_redundant_cut
		!quiet && @warn "The Redundant-Cut is only used when faithful2furini2016" *
			" is enabled. As flag faithful2furini2016 is not enabled, flag" *
			" no-redundant-cut has no effect."
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
	#plis = Dict{Tuple{S, S}, P}()
	plis = zeros(P, L, W)
	#plis[(L, W)] = one(P)
	plis[L, W] = one(P)
	# next: plates already indexed but not yet processed, starts with (L, W, 1).
	# Storing the index as the third value is not necessary (as it could be
	# queried from plis) but this is probably more efficient this way.
	next = Vector{Tuple{S, S, P}}()
	next_idx = one(P)
	#sizehint!(next, max_num_plates)
	push!(next, (L, W, one(P)))
	if !no_redundant_cut
		# If the preprocessing is faithful2furini2016 and Redundant-Cut is
		# enabled, then four auxiliary trim cut flag vectors are necessary.
		sfhv = (Vector{Int8}(), Vector{Int8}(), Vector{Int8}(), Vector{Int8}())
		sh, sv, fh, fv = sfhv
		push!.(sfhv, (1, 1, 1, 1))
	end
	# n: The amount of plates (the index of the highest plate type).
	n = one(P) # there is already the original plate
	# Memoized discretizations. The discretized lengths for every plate width.
	# The discretized widths for every plate length.
	dls = [Vector{S}() for _ = 1:W]
	dws = [Vector{S}() for _ = 1:L]
	# There are two discretizations for each dimension: dl1/dw1 and dl2/dw2.
	# They are only different if the Cut-Position reduction is being used. In
	# this case, dl1/dw1, which is used for cutting the plate, may be a
	# restricted cut set (i.e., it does not have piece size combinations, just
	# single piece sizes). The dl2/dw2 is always the complete discretization,
	# and may be needed (even if Cut-Position is used) for reducing the size of
	# the second child to the last discretized position (if round2disc is
	# enabled).
	dl1 = dl2 = dls[W] = discretize(
		d, l, w, L, W; ignore_W = ignore_2th_dim, ignore_d = ignore_d
	)
	#@show length(dl1)
	dw1 = dw2 = dws[L] = discretize(
		d, w, l, W, L; ignore_W = ignore_2th_dim, ignore_d = ignore_d
	)
	#fpr_count = 0
	#@show length(dw1)
	while next_idx <= length(next)
		#@show next_idx
		#@show length(next)
		@assert (no_redundant_cut || [length(next)] == unique(length.(sfhv)))

		pll, plw, pli = next[next_idx] # PLate Length, Width, and Index
		if !no_redundant_cut
			sh_j, sv_j, fh_j, fv_j = getindex.(sfhv, next_idx)
		end
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
			isempty(dls[plw]) && (dls[plw] = discretize(
				d, l, w, L, plw; ignore_d = ignore_d
			))
			dl1 = dl2 = dls[plw]

			isempty(dws[pll]) && (dws[pll] = discretize(
				d, w, l, W, pll; ignore_d = ignore_d
			))
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

			# If the preprocessing is faithful2furini2016, the plates are cut until
			# they have the same size as pieces and, consequently, there exist the
			# concept of trim cut (i.e., a cut in which the second child is waste).
			trim_cut = faithful2furini2016 && !fits_at_least_one(sllw, pll, plw - x)
			# The trim_cut flag is used below outside of 'if faithful2furini2016'
			# because a true value implicates that 'faithful2furini2016 == true'.
			#if iszero(get(plis, (pll, x), 0)) # If the first child does not yet exist.
			if iszero(plis[pll, x]) # If the first child does not yet exist.
				push!(next, (pll, x, n += 1)) # Create the first child.
				if !no_redundant_cut
					if trim_cut
						# From Furini2016 supplement: "if a NEW plate j_1 ∈ J is obtained
						# from j through a trim cut with orientation v:"
						fchild_sfhv = (0, -1, sh_j, sv_j)
					else
						# If a plate (NEW or existing) j_1 is obtained from j without a trim
						# cut: set all flags of j 1 to 1.
						fchild_sfhv = (1, 1, 1, 1)
					end
					push!.(sfhv, fchild_sfhv)
				end
				#plis[(pll, x)] = n # Mark plate existence.
				plis[pll, x] = n # Mark plate existence.
			elseif !no_redundant_cut
				# If plate already exists, and we are implementing Redundant-Cut
				if trim_cut
					# From Furini2016 supplement: "if an existing plate j_1 ∈ J is
					# obtained from j through a trim cut:"
					#sh_j > -1 && (fh[plis[(pll, x)]] = 1)
					#sv_j > -1 && (fv[plis[(pll, x)]] = 1)
					sh_j > -1 && (fh[plis[pll, x]] = 1)
					sv_j > -1 && (fv[plis[pll, x]] = 1)
				else
					# From Furini2016 supplement: "if a plate (new or EXISTING) j_1 is
					# obtained from j without a trim cut: set all flags of j_1 to 1."
					#setindex!.(sfhv, (1, 1, 1, 1), plis[(pll, x)])
					setindex!.(sfhv, (1, 1, 1, 1), plis[pll, x])
				end
			end
			# If the second child size is rounded down to a discretized point:
			if round2disc && !trim_cut
				# Takes the index of the rounded down discretized point, if it is zero,
				# it does not exist. Such may happen if faithful2furini2016 is enabled,
				# in the case the second child is waste (smaller than the first
				# discretized point).
				dw_sc_ix = searchsortedlast(dw2, plw - x)
				@assert faithful2furini2016 || !iszero(dw_sc_ix)
				if iszero(dw_sc_ix)
					dw_sc = zero(S)
				else
					dw_sc = dw2[dw_sc_ix]
				end
			else # If the second child size is not rounded.
				dw_sc = plw - x
			end

			# If the second child is not waste (!trim_cut), nor was already
			# generated, then it is a new plate obtained without a trim cut.
			#if !trim_cut && iszero(get(plis, (pll, dw_sc), 0))
			if !trim_cut && iszero(plis[pll, dw_sc])
				# Assert meaning: if the second child is not waste, then it must
				# have an associated size.
				@assert !iszero(dw_sc)
				# Save the plate to be processed later, and mark its existence.
				push!(next, (pll, dw_sc, n += 1))
				#plis[(pll, dw_sc)] = n
				plis[pll, dw_sc] = n
				# From Furini2016 supplement: "if a plate (NEW or existing) j 1 is
				# obtained from j without a trim cut: set all flags of j_1 to 1."
				!no_redundant_cut && push!.(sfhv, (1, 1, 1, 1))
			end
			# Add the cut to the cut list. Check if the second child plate is waste
			# before trying to get its plate index.
			#push!(vnnn, (pli, plis[(pll, x)], trim_cut ? 0 : plis[(pll, dw_sc)]))
			push!(vnnn, (pli, plis[pll, x], trim_cut ? 0 : plis[pll, dw_sc]))
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

			# If the preprocessing is faithful2furini2016, the plates are cut until
			# they have the same size as pieces and, consequently, there exist the
			# concept of trim cut (i.e., a cut in which the second child is waste).
			trim_cut = faithful2furini2016 && !fits_at_least_one(sllw, pll - y, plw)
			# The trim_cut flag is used below outside of 'if faithful2furini2016'
			# because a true value implicates that 'faithful2furini2016 == true'.
			#if iszero(get(plis, (y, plw), 0)) # If the first child does not yet exist.
			if iszero(plis[y, plw]) # If the first child does not yet exist.
				push!(next, (y, plw, n += 1)) # Create the first child.
				if !no_redundant_cut
					if trim_cut
						# From Furini2016 supplement: "if a NEW plate j_1 ∈ J is obtained
						# from j through a trim cut with orientation v:"
						fchild_sfhv = (-1, 0, sh_j, sv_j)
					else
						# If a plate (NEW or existing) j_1 is obtained from j without a trim
						# cut: set all flags of j 1 to 1.
						fchild_sfhv = (1, 1, 1, 1)
					end
					push!.(sfhv, fchild_sfhv)
				end
				#plis[(y, plw)] = n # Mark plate existence.
				plis[y, plw] = n # Mark plate existence.
			elseif !no_redundant_cut
				# If plate already exists, and we are implementing Redundant-Cut
				if trim_cut
					# From Furini2016 supplement: "if an existing plate j_1 ∈ J is
					# obtained from j through a trim cut:"
					#sh_j > -1 && (fh[plis[(y, plw)]] = 1)
					#sv_j > -1 && (fv[plis[(y, plw)]] = 1)
					sh_j > -1 && (fh[plis[y, plw]] = 1)
					sv_j > -1 && (fv[plis[y, plw]] = 1)
				else
					# From Furini2016 supplement: "if a plate (new or EXISTING) j_1 is
					# obtained from j without a trim cut: set all flags of j_1 to 1."
					#setindex!.(sfhv, (1, 1, 1, 1), plis[(y, plw)])
					setindex!.(sfhv, (1, 1, 1, 1), plis[y, plw])
				end
			end
			# If the second child size is rounded down to a discretized point:
			if round2disc && !trim_cut
				# Takes the index of the rounded down discretized point, if it is zero,
				# it does not exist. Such may happen if faithful2furini2016 is enabled,
				# in the case the second child is waste (smaller than the first
				# discretized point).
				dl_sc_ix = searchsortedlast(dl2, pll - y)
				@assert faithful2furini2016 || !iszero(dl_sc_ix)
				if iszero(dl_sc_ix)
					dl_sc = zero(S)
				else
					dl_sc = dl2[dl_sc_ix]
				end
			else # If the second child size is not rounded.
				dl_sc = pll - y
			end
			# If the second child is not waste (!trim_cut), nor was already
			# generated, then it is a new plate obtained without a trim cut.
			#if !trim_cut && iszero(get(plis, (dl_sc, plw), 0))
			if !trim_cut && iszero(plis[dl_sc, plw])
				# Assert meaning: if the second child is not waste, then it must
				# have an associated size.
				@assert !iszero(dl_sc)
				# Save the plate to be processed later, and mark its existence.
				push!(next, (dl_sc, plw, n += 1))
				#plis[(dl_sc, plw)] = n
				plis[dl_sc, plw] = n
				# From Furini2016 supplement: "if a plate (NEW or existing) j 1 is
				# obtained from j without a trim cut: set all flags of j_1 to 1."
				!no_redundant_cut && push!.(sfhv, (1, 1, 1, 1))
			end
			# Add the cut to the cut list. Check if the second child plate is waste
			# before trying to get its plate index.
			#push!(hnnn, (pli, plis[(y, plw)], trim_cut ? 0 : plis[(dl_sc, plw)]))
			push!(hnnn, (pli, plis[y, plw], trim_cut ? 0 : plis[dl_sc, plw]))
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
			@assert !iszero(plis[l[pii], w[pii]])
			push!(np, (plis[l[pii], w[pii]], pii))
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

	first_vertical_cut_idx = length(hnnn) + 1
	return ByproductPPG2KP(
		append!(hnnn, vnnn), first_vertical_cut_idx, np, pli_lwb,
		deepcopy(d), deepcopy(l), deepcopy(w), deepcopy(L), deepcopy(W)
	)
end

end # module

