module Enumeration

# TODO: Consider removing flag ignore_2th_dim. Reasoning: nobody does that,
# even Furini only consider the pieces that fit a plate before discretizing it,
# even if ignoring the 2th dimension would save a lot of effort in the
# preprocessing nobody believes would be a valid trade-off to do so (the
# preprocessing phase participation on the total time is negligible), and it
# makes the code far more complex (because I have to have code that respects
# and ignores a basic assumption).

export SortedLinkedLW
export becker2019_discretize, gen_cuts, gen_cuts_sb
export should_extract_piece_from_plate

# Structure for keeping the piece lengths and widths both in the original order
# and sorted order, and allow to access the piece index or the width (length)
# counterpart when iterating the sorted lengths (widths).
struct SortedLinkedLW{D, S}
	l :: Vector{S}
	w :: Vector{S}
	sl :: Vector{S}
	sw :: Vector{S}
	sli2pii :: Vector{D} # sorted length index to piece index
	swi2pii :: Vector{D} # sorted width index to piece index
	pii2sli :: Vector{D} # piece index to sorted length index
	pii2swi :: Vector{D} # piece index to sorted width index
end

# IMPORTANT: this constructor does not copy l and w, so any changes to l or w
# will reflect in it, silently and completely invalidating the structure.
# The first parameter is the type of integer to be used in the indexes.
# This constructor basically creates a sorted copy of both vectors, then
# create the four inverse indexes, and return them all in a single structure.
# This constructor is not implemente in the most efficient way but it is
# called one single time for each instance, so, irrelevant.
function SortedLinkedLW(::Type{D}, l :: Vector{S}, w :: Vector{S}) where {D, S}
	@assert length(l) == length(w)
	n = length(l)
	sl = sort(l)
	sw = sort(w)
	sli2pii = sort!(collect(1:n), by = pii -> l[pii])
	swi2pii = sort!(collect(1:n), by = pii -> w[pii])
	pii2sli = Vector{D}(undef, n)
	pii2swi = Vector{D}(undef, n)
	for si = 1:n
		pii2sli[sli2pii[si]] = si
		pii2swi[swi2pii[si]] = si
	end
	SortedLinkedLW(l, w, sl, sw, sli2pii, swi2pii, pii2sli, pii2swi)
end

# TODO: break this method into two (like it was done for the flow model)
#       one for reducing the input, other to discretize only over one dimension
# Intelligent discretization method that: (1) takes demand into account;
# (2) takes the dimension that is not being discretized into account. It can
# also be used to discover which lengths are only shared by single pieces
# and not linear combinations of two or more pieces.
# d: piece demand (d)
# l: piece dimensions for dim that is being discretized (may be `l` or `w`)
# w: piece dimensions for dim that is NOT being discretized (may be `l` or `w`)
# L: plate dimension for dim that is being discretized (may be `L` or `W)`
# W: plate dimension for dim that is NOT being discretized (may be `L` or `W`)
# mark_single_piece: instead of discretizing the plate dimension, give all
#   the piece dimensions (for the dimension being discretized) that are shared
#   only by single pieces and not by two or more pieces (even of the same type)
#   combinations.
function becker2019_discretize(
	d :: Vector{D}, l :: Vector{S}, w :: Vector{S}, L :: S, W :: S;
	only_single_pieces = false, ignore_W = false, ignore_d = false
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

# Check if such piece should be extracted from a plate. This method may be used
# from both symmetry-breaking and non-symmetry-breaking contexts, because
# the default value of parameter symm is the allow-any-cut-direction value.
# A piece should be extracted from a plate if:
# (1) The plate allows the piece size.
# (2) The plate will not be cut in a way that leaves a second child large
#     enough to allow the piece size (i.e., if the piece extraction may be
#     postponed to be done in a second child of the current plate, it will be
#     postponed).
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

# TODO: internal change. The symmetry tag dimension of plis should be the
# first one, given how columnwise memory layout works.
function gen_cuts_sb(
	::Type{P}, d :: Vector{D}, sllw :: SortedLinkedLW{D, S}, L :: S, W :: S;
	ignore_2th_dim = false, ignore_d = false, round2disc = true
) where {D, S, P}
	(ignore_d || ignore_2th_dim || round2disc) && @error "ignore_2th_dimm, ignore_d, and round2disc are not yet implemented for gen_cuts_sb, also, first improve the performance by memoizing the discretizations"
	l = sllw.l
	w = sllw.w
	@assert length(d) == length(l)
	@assert length(d) == length(w)
	only_single_l = becker2019_discretize(
		d, l, w, L, W; only_single_pieces = true, ignore_W = ignore_2th_dim,
		ignore_d = ignore_d
	)
	only_single_w = becker2019_discretize(
		d, w, l, W, L; only_single_pieces = true, ignore_W = ignore_2th_dim,
		ignore_d = ignore_d
	)
	max_piece_type = convert(D, length(l))
	max_num_plates = convert(P, L) * convert(P, W) * 3
	min_pil = minimum(l)
	min_piw = minimum(w)
	# If (n1, n2, n3) is in hcut (vcuts), then n1 may be partitioned in n2 and n3
	# by the means of a horizontal (vertical) cut. Not every possible partition
	# is present, just the ones which each child plate can fit at least one
	# piece, and the ones where one child plate has the same dimension as a piece
	# and the other is the dummy waste plate.
	hcuts = Vector{NTuple{3, P}}()
	vcuts = Vector{NTuple{3, P}}()
	# If the plate has the exact same size that a piece, then it is added to
	# the vector at pii2plis[piece_type].
	pii2plis = [Vector{P}() for _ = 1:max_piece_type]
	# If the piece has the exact same size that a plate, then it is added to
	# the vector at pli2piis[plate_type].
	pli2piis = Vector{D}[]
	# The list of plates attributes: plate length, plate width, plate symmetry,
	# and plate bound. The plate index is the same as the index in pli_lwsb.
	pli2lwsb = Vector{Tuple{S, S, UInt8, P}}()
	# plis: matrix of the plate dimensions in which zero means "never seen that
	# plate before" and nonzero means "this nonzero number is the plate index".
	plis = zeros(P, L, W, 3)
	plis[L, W, 3] = one(P)
	# next: plates already indexed but not yet processed, starts with (L, W,
	# 3, 1). The size of the original plate, the code for being able to cut
	# in any direction (same as third dimension of plis), and the plate index.
	# Storing the index as the third value is not necessary (as it could be
	# queried from plis) but this is probably more efficient this way.
	next = Vector{Tuple{S, S, UInt8, P}}()
	next_idx = one(P)
	#sizehint!(next, max_num_plates)
	push!(next, (L, W, 3, one(P)))
	# n: The amount of plates (the index of the highest plate type).
	n = one(P) # there is already the original plate
	piece_cuts_avoided_by_symmb = 0
	while next_idx <= length(next)
		# PLate Length, Width, Symmetry, and Index
		pll, plw, pls :: UInt8, pli = next[next_idx]
		next_idx += 1
		plb = (L ÷ pll) * (W ÷ plw) # PLate Bound
		# It is not necessary to store the plate id in pli2lwsb because they are
		# added in order (a queue is used), so the array index is the plate index.
		push!(pli2lwsb, (pll, plw, pls, plb))
		# the pli2piis vector grows with the number of plates processed
		push!(pli2piis, D[])
		for pii in 1:max_piece_type # pii: PIece Index
			pil, piw = l[pii], w[pii] # PIece Length and Width (homophones, I know)
			if pil == pll && piw == plw
				# If the current plate is exactly the size of a plate, add it to
				# the list of plate codes that correspond to a piece.
				push!(pii2plis[pii], pli)
				push!(pli2piis[pli], pii)
			elseif should_extract_piece_from_plate(pii, pll, plw, sllw, pls)
				# If there is no way to cut the current piece from the current plate
				# except if we cut after the half of the plate... (i.e., the most
				# unbalanced cut possible, that is, making a cut with the minimal
				# length or width among the pieces, removes the possibility of getting
				# such piece from the larger child plate)
				if pil < pll && (pls == 1 || pls == 3)
					# If the current plate is allowed to be cut horizontally (i.e.,
					# reducing the length) and we need to reduce the length to make the
					# plate into the piece, then we do so.
					if iszero(plis[pil, plw, 2])
						push!(next, (pil, plw, 2, n += 1))
						plis[pil, plw, 2] = n
					end
					push!(hcuts, (pli, plis[pil, plw, 2], 0))
				elseif piw < plw && (pls == 2 || pls == 3)
					# If the current plate is allowed to be cut vertically (i.e.,
					# reducing the width) and we need to reduce the width to make the
					# plate into the piece, then we do so.
					if iszero(plis[pll, piw, 1])
						push!(next, (pll, piw, 1, n += 1))
						plis[pll, piw, 1] = n
					end
					push!(vcuts, (pli, plis[pll, piw, 1], 0))
				else
					piece_cuts_avoided_by_symmb += 1
				end
			end
		end
		# TODO: after is working, refactor the code below to use a single loop
		# using a vector saved in a variable, or an inner function
		sorted_in(a, v) = searchsortedfirst(a, v) <= length(a)
		if pls == 1 && sorted_in(only_single_w, plw)
			indexes = filter(i -> w[i] == plw && l[i] < pll, collect(1:max_piece_type))
			#@assert !isempty(indexes)
			for y in l[indexes]
				if iszero(plis[y, plw, 2])
					push!(next, (y, plw, 2, n += 1))
					plis[y, plw, 2] = n
				end
				if iszero(plis[pll - y, plw, 3])
					push!(next, (pll - y, plw, 3, n += 1))
					plis[pll - y, plw, 3] = n
				end
				push!(hcuts, (pli, plis[y, plw, 2], plis[pll - y, plw, 3]))
			end
		elseif pls == 1 || pls == 3
			for y in becker2019_discretize(
				d, l, w, pll ÷ 2, plw; ignore_2th_dim = ignore_2th_dim
			)
				#y > pll ÷ 2 && break
				@assert plw >= min_piw
				@assert y >= min_pil
				@assert pll - y >= min_pil
				if iszero(plis[y, plw, 2])
					push!(next, (y, plw, 2, n += 1))
					plis[y, plw, 2] = n
				end
				if iszero(plis[pll - y, plw, 3])
					push!(next, (pll - y, plw, 3, n += 1))
					plis[pll - y, plw, 3] = n
				end
				push!(hcuts, (pli, plis[y, plw, 2], plis[pll - y, plw, 3]))
			end
		end
		if pls == 2 && sorted_in(only_single_l, pll)
			indexes = filter(i -> l[i] == pll && w[i] < plw, collect(1:max_piece_type))
			#@assert !isempty(indexes)
			for x in w[indexes]
				if iszero(plis[pll, x, 1])
					push!(next, (pll, x, 1, n += 1))
					plis[pll, x, 1] = n
				end
				if iszero(plis[pll, plw - x, 3])
					push!(next, (pll, plw - x, 3, n += 1))
					plis[pll, plw - x, 3] = n
				end
				push!(vcuts, (pli, plis[pll, x, 1], plis[pll, plw - x, 3]))
			end
		elseif pls == 2 || pls == 3
			for x in becker2019_discretize(
				d, w, l, plw ÷ 2, pll; ignore_2th_dim = ignore_2th_dim
			)
				#x > plw ÷ 2 && break
				@assert pll >= min_pil
				@assert x >= min_piw
				@assert plw - x >= min_piw
				if iszero(plis[pll, x, 1])
					push!(next, (pll, x, 1, n += 1))
					plis[pll, x, 1] = n
				end
				if iszero(plis[pll, plw - x, 3])
					push!(next, (pll, plw - x, 3, n += 1))
					plis[pll, plw - x, 3] = n
				end
				push!(vcuts, (pli, plis[pll, x, 1], plis[pll, plw - x, 3]))
			end
		end
	end
	#@show piece_cuts_avoided_by_symmb

	# For each possible subplatex size, same_size_plis groups the plate indexes
	# that refer to the plates with the same size, and distinct symmetry tags.
	same_size_plis = Vector{Vector{P}}()
	# The length of same_size_plis may be as low as one third of length(pli2lwsb)
	# and as high as length(pli2lwsb).
	sizehint!(same_size_plis, length(pli2lwsb))
	for j = 1:W, i = 1:L # column-major for performance
		temp = Vector{P}()
		!iszero(plis[i, j, 1]) && push!(temp, plis[i, j, 1])
		!iszero(plis[i, j, 2]) && push!(temp, plis[i, j, 2])
		!iszero(plis[i, j, 3]) && push!(temp, plis[i, j, 3])
		!iszero(length(temp)) && (push!(same_size_plis, temp); temp = Vector{P}())
	end
	@assert last(same_size_plis) == P[1]
	pop!(same_size_plis)

	return pli2lwsb, hcuts, vcuts, pii2plis, pli2piis, same_size_plis
end

# Take pieces and a plate, return an upper bound on the number of pieces
# that fit inside the plate. If cutoff becomes greater than limit, return
# immediately (does not provide a real upper bound, you just know the real
# upper bound is greater than cutoff). If ignore_2th_dim is true, then it
# does not use the information of w and W.
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
				ub > cutoff && return ub
			end
		end
	end
	return ub
end

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
		d, l, w, a, L, W, A, 5; ignore_2th_dim = ignore_2th_dim
	)
	return ub <= 5
end

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

function fits_at_least_one(
	sllw :: SortedLinkedLW{D, S}, L :: S, W :: S
) :: Bool where {D, S}
	for (i, v) in enumerate(sllw.sl)
		v > L && return false
		sllw.w[sllw.sli2pii[i]] <= W && return true
	end

	return false
end

function filter_symm_pos(
	disc :: Vector{S},
	dim :: S,
	disc_copy = copy(disc) :: Vector{S}
) :: Vector{S} where {S}
	half_dim = dim ÷ 2
	for xy in disc
		# If the plate has
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

function gen_cuts(
	::Type{P}, d :: Vector{D}, sllw :: SortedLinkedLW{D, S}, L :: S, W :: S;
	ignore_2th_dim = false, ignore_d = false, round2disc = true,
	faithful2furini2016 = false,
	no_redundant_cut = false, no_cut_position = false,
	no_furini_symmbreak = false
) where {D, S, P}
	faithful2furini2016 && round2disc && @warn(
		"Enabling both faithful2furini2016 and round2disc is allowed, but you are not being entirely faithful to Furini2016 if you do so."
	)
	!faithful2furini2016 && no_redundant_cut && @warn(
		"The Redundant-Cut is only used when faithful2furini2016 is enabled. As flag faithful2furini2016 is not enabled, flag no-redundant-cut has no effect."
	)
	!faithful2furini2016 && (no_redundant_cut = true)
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
	dl1 = dl2 = dls[W] = becker2019_discretize(
		d, l, w, L, W; ignore_W = ignore_2th_dim, ignore_d = ignore_d
	)
	#@show length(dl1)
	dw1 = dw2 = dws[L] = becker2019_discretize(
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
			isempty(dls[plw]) && (dls[plw] = becker2019_discretize(
				d, l, w, L, plw; ignore_d = ignore_d
			))
			dl1 = dl2 = dls[plw]

			isempty(dws[pll]) && (dws[pll] = becker2019_discretize(
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
						fchild_sfhv = (-1, 0, sh_j, sv_j)
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
						fchild_sfhv = (0, -1, sh_j, sv_j)
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
			filter_redundant_cuts!(hnnn, sh, fh)
			filter_redundant_cuts!(vnnn, sv, fv)
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

	return pli_lwb, hnnn, vnnn, np
end

end # module

