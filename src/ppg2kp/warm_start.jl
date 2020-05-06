# Returns the length (width) of the piece of smallest length (width) that fits
# the plate defined by LxW. Errors if there is no piece that fits the plate.
function _min_l_of_fit_piece(
	sllw :: SortedLinkedLW{D, S}, L :: S, W :: S
) :: S where {D, S}
	for (i, l) in enumerate(sllw.sl)
		l > L && break
		sllw.w[sllw.sli2pii[i]] <= W && return l
	end
	@error "No piece fits the plate"
end
function _min_w_of_fit_piece(
	sllw :: SortedLinkedLW{D, S}, L :: S, W :: S
) :: S where {D, S}
	for (i, w) in enumerate(sllw.sw)
		w > W && break
		sllw.l[sllw.swi2pii[i]] <= L && return w
	end
	@error "No piece fits the plate"
end

# Returns a vector with length equal to highest_plate, and in which position
# `p` gives the range of `cuts` that have `p` as their first field (i.e., the
# parent plate). This range is shifted by offset.
# NOTE: `highest_plate >= maximum(getindex.(cuts, 1))` or undefined behavior.
# NOTE: if some plate `p` (1 < p <= highest_plate) is not the first field of
# any element in `cuts`, then the range of `p` in the returned vector will be
# an empty range and, consequently, the returned array will probably fail a
# `issorted` check, because no empty range is considered greater than a
# non-empty range;
function _build_pp_ranges(
	cuts :: AbstractVector{NTuple{3, P}},
	highest_plate :: P,
	offset :: P = zero(P)
) :: Vector{UnitRange{P}} where {P}
	# If `cuts` is empty, return a vector of empty ranges of the expected size.
	empty_r = one(P):zero(P)
	isempty(cuts) && return fill(empty_r, highest_plate)
	# If there is a cut, there should be at least one plate.
	@assert highest_plate >= one(P)
	ranges = Vector{UnitRange{P}}(undef, highest_plate)
	curr_pp = first(cuts)[PARENT]
	# The original plate (index one) is always present.
	@assert isone(curr_pp)
	first_idx_curr_pp = one(P)
	# The first field of `cuts` elements (i.e., PARENT) never decreases, and when
	# it increases, the range in which they kept the same value is saved.
	for (i, cut) in pairs(cuts)
		pp = cut[PARENT]
		if pp < curr_pp
			@error "The parameter `cuts` must be sorted by its first field."
		elseif pp > curr_pp
			ranges[curr_pp] = (offset + first_idx_curr_pp):(offset + i - 1)
			# If pp - curr_pp > 1 (i.e., pp values were skipped) save empty
			# ranges for such values.
			for j = (curr_pp+1):(pp-1); ranges[j] = empty_r; end
			# Update our sentinels.
			curr_pp = pp
			first_idx_curr_pp = i
		end
		# If `curr_pp` is equal to `highest_plate` we assume `cuts` is
		# non-decreasing and break early. The `ranges[highest_plate]` is
		# saved after the loop.
		curr_pp == highest_plate && break
	end
	# The last `curr_pp` has no change of value to trigger its writting.
	ranges[curr_pp] = (offset + first_idx_curr_pp):(offset + length(cuts))
	return ranges
end

function _get_LW(
	bp :: ByproductPPG2KP{D, S, P},
	plate_idx :: P
) :: Tuple{S, S} where {D, S, P}
	return bp.pli_lwb[plate_idx][LENGTH], bp.pli_lwb[plate_idx][WIDTH]
end

# Find a cut in `bp.cuts[cuts_by_pp[pp]]` with a FIRST_CHILD of size `cut_size`
# in the dimension `dim`, `push!` it to `cuts_done`, and return its
# `FIRST_CHILD` and `SECOND_CHILD` (a convenience, as we could use the last
# element of `cuts_done` to index `bp.cuts` and get them).
function _make_cut!(
	cuts_done :: Vector{P},
	pp :: P,
	bp :: ByproductPPG2KP{D, S, P},
	cuts_by_pp :: Vector{UnitRange{P}},
	cut_size :: S,
	dim :: Dimension
) where {D, S, P}
	idx_in_view = findfirst(
		c -> bp.pli_lwb[c[FIRST_CHILD]][dim] == cut_size,
		view(bp.cuts, cuts_by_pp[pp])
	) :: Int
	cut_idx = cuts_by_pp[pp].start + convert(P, idx_in_view) - one(P)
	push!(cuts_done, cut_idx)
	@assert bp.cuts[cut_idx][PARENT] == pp
	return bp.cuts[cut_idx][FIRST_CHILD], bp.cuts[cut_idx][SECOND_CHILD]
end

# For `qt_trims` times, find a cut that reduces the current plate by
# `trim_size` (the method errors if such cut does not exist) and update
# the current plate, return the remaining plate after all those trim cuts,
# saves the trim cuts done to cuts_done.
function _successive_trims!(
	cuts_done :: Vector{P},
	plate_idx :: P,
	bp :: ByproductPPG2KP{D, S, P},
	cuts_by_pp :: Vector{UnitRange{P}},
	qt_trims :: P,
	trim_size :: S,
	dim :: Dimension
) :: P where {D, S, P}
	pp = plate_idx
	for _ = 1:qt_trims
		# We throw away the info about which was the first child (i.e., the trim
		# plate) and update the parent plate to be the second child.
		_, pp = _make_cut!(cuts_done, pp, bp, cuts_by_pp, trim_size, dim)
		@assert !iszero(pp)
	end
	return pp
end

# Trim the dimension `dim` of plate `pp` to size `s` (if `faithful2furini`) or
# close enough (if `!faithful2furini`, and by 'close enough' we mean that no
# piece would fit in the trim of the 'close enough' plate to exactly size `s`).
# The cuts needed are added to `cuts_done`. For a plate index `pli` the
# `cuts_by_pp[pli]` has the range of cuts (that divide `dim`) with `pli` as
# the PARENT plate of the cut.
function _safe_trim_dim!(
	cuts_done :: Vector{P},
	pp :: P,
	s :: S,
	dim :: Dimension,
	cuts_by_pp :: Vector{UnitRange{P}},
	sllw :: SortedLinkedLW{D, S},
	bp :: ByproductPPG2KP{D, S, P},
	faithful2furini :: Bool
) where {D, S, P}
	# rS: size of `dim` in plate `pp`
	rS = bp.pli_lwb[pp][dim]
	@assert rS >= s
	if s != rS # Check if trimming is actually necessary.
		# If the cut is in the first half of the plate it is guaranteed to exist.
		if s <= (rS+1)/2
			# In the case of a 'clean' trim cut, the trim is the SECOND_CHILD that
			# is thrown away, and the PARENT is updated to the FIRST_CHILD.
			pp, _ = _make_cut!(cuts_done, pp, bp, cuts_by_pp, s, dim)
		else # If a trim cut may be needed.
			L, W = _get_LW(bp, pp)
			if dim == LENGTH
				min_s = _min_l_of_fit_piece(sllw, L, W)
			else
				min_s = _min_w_of_fit_piece(sllw, L, W)
			end
			# If the difference between the plate and the piece allows for copies of
			# the smallest-dim-size piece-that-fits, reduce the difference by
			# successive trim cuts.
			if rS - s >= min_s
				qt_trims = convert(P, div(rS - s, min_s))
				pp = _successive_trims!(
					cuts_done, pp, bp, cuts_by_pp, qt_trims, min_s, dim
				)
				rS = bp.pli_lwb[pp][dim]
				# If !faithful2furini and the other dim has a difference smaller than
				# the smallest-other-dim-size piece-that-fits, then it is guaranteed
				# that the extraction will be in bp.np.
			end
			# Otherwise, if faithful2furini, the match of piece and plate need
			# to be exact, so we may need a final trim cut. It was done here,
			# after the successive trim cuts, because here it is guaranteed to
			# exist. If it was done before, it could be the cut did not exist
			# because the symmetry removal used by Furini2016 (i.e, if there was a
			# cut of length exactly L - l, in the first half of the plate, then a
			# cut of length l, in the second half of the plate, would not exist).
			if faithful2furini && s != rS
				pp, _ = _make_cut!(cuts_done, pp, bp, cuts_by_pp, s, dim)
			end
		end
	end
	return pp
end

function _sell_piece!(
	pieces_sold :: Vector{P},
	cuts_done :: Vector{P},
	plate_idx :: P,
	piece_idx :: D,
	sllw :: SortedLinkedLW{D, S},
	bp :: ByproductPPG2KP{D, S, P},
	hcuts_by_pp :: Vector{UnitRange{P}},
	vcuts_by_pp :: Vector{UnitRange{P}},
	faithful2furini :: Bool
) :: Nothing where {D, S, P}
	# The piece-to-be-sold dimensions.
	piece_l, piece_w = bp.l[piece_idx], bp.w[piece_idx]
	# pp: the current parent plate that is updated as we trim the given plate
	pp = plate_idx
	pp = _safe_trim_dim!(
		cuts_done, pp, piece_l, LENGTH, hcuts_by_pp, sllw, bp, faithful2furini
	)
	pp = _safe_trim_dim!(
		cuts_done, pp, piece_w, WIDTH, vcuts_by_pp, sllw, bp, faithful2furini
	)
	# Now a extraction from the current plate to piece_idx MUST exist.
	np_idx = findfirst(==((pp, piece_idx)), bp.np) :: Int
	push!(pieces_sold, convert(P, np_idx))
	return
end

# Given a 2-staged `pattern` and the preprocessed information of a PPG2KP model
# (i.e., `bp` the byproduct), a vector of indexes of bp.cuts and a vector
# of indexes of bp.np that together represent a valid solution of the PPG2KP
# model representing the solution given by the 2-staged cut pattern. This is
# used to be able to warm start the model from an heuristic solution.
"""
    cuts_and_extractions_from_2_staged_solution(pattern, bp, faithful2furini)

Given a 2-staged `pattern` and a PPG2KP model `bp` (byproduct), returns the
cuts (`bp.cuts` indexes) and extractions (`bp.np` indexes) that make up a
valid solution representing such `pattern`.

The `pattern` parameter has many restrictions:

1) It is 2-staged, what is enforced by being a vector of vectors.
2) Each inner vector is a vertical stripe with the width of the first piece
   inside such inner vector. No horizontal stripes allowed. No piece inside
   an inner vector has width larger than the first piece of the inner vector.
3) The outer vector is sorted by non-increasing stripe width.
4) No inner vectors should be empty. Every number is a valid piece type index.
"""
function cuts_and_extractions_from_2_staged_solution(
	pattern :: Vector{Vector{D}},
	bp :: ByproductPPG2KP{D, S, P},
	faithful2furini :: Bool
) :: Tuple{Vector{P}, Vector{P}} where {D, S, P}
	cuts_done, pieces_sold = P[], P[]
	isempty(pattern) && return cuts_done, pieces_sold
	@assert all((!isempty).(pattern)) # No stripe should be empty.
	# @assert explanation: the first stripe is the one of largest width.
	@assert all(fp -> bp.w[fp] <= bp.w[pattern[1][1]], getindex.(pattern, 1))
	# @assert explanation: each of the inner vectors has no piece with
	# width greater than their width (i.e., the width of the first element).
	@assert all(a -> all(p -> bp.w[p] <= bp.w[a[1]], a), pattern)
	sllw = SortedLinkedLW(D, bp.l, bp.w)
	highest_plate = length(bp.pli_lwb)
	# Compute cut ranges that share the same orientation (vertical or horizontal)
	# and the same parent plate (pp). This is O(n) and greatly reduces the effort
	# spent finding the desired cuts in the rest of the algorithm.
	hcuts_by_pp = _build_pp_ranges(
		(@view bp.cuts[one(P):(bp.first_vertical_cut_idx - one(P))]), highest_plate
	)
	num_cuts = length(bp.cuts)
	vcuts_by_pp = _build_pp_ranges(
		(@view bp.cuts[bp.first_vertical_cut_idx:num_cuts]), highest_plate,
		bp.first_vertical_cut_idx - one(P)
	)
	rp = one(P) # The initial `r`emaining `p`late is the original plate.
	num_stripes = length(pattern)
	# stripe_plates: store the plate for each stripe, used after but not updated.
	stripe_plates = Vector{P}(undef, num_stripes)
	# The order we make the cuts is not specially relevant, except that we only
	# have the guarantee that the cuts really exist if they are in the first
	# half of the plate, and this happens naturally if we left the largest
	# width stripe for last.
	for vstripe_idx in 2:num_stripes
		# Cut a vertical (i.e., keep the length divides the width) stripe.
		stripe_plates[vstripe_idx], rp = _make_cut!(
			cuts_done, rp, bp, vcuts_by_pp, bp.w[pattern[vstripe_idx][1]], WIDTH
		)
	end
	# If it exists a clean cut for the last stripe do it.
	largest_w = bp.w[pattern[1][1]]
	if largest_w <= (bp.pli_lwb[rp][WIDTH]+1)/2
		stripe_plates[1], _ = _make_cut!(
			cuts_done, rp, bp, vcuts_by_pp, largest_w, WIDTH
		)
	else # Otherwise, just use the remaining plate as last stripe.
		stripe_plates[1] = rp
	end
	# Now all the vertical stripes are cut. The stripes need to be cut
	# horizontally creating a new plate for each piece inside the stripe,
	# and then these plates may need a vertical trim (because they are smaller
	# than the stripe width) before they are finally sold as pieces.
	rp = zero(P) # From there on, only the `stripe_rp` is relevant.
	for vstripe_idx in 1:num_stripes # the outer order is not relevant
		stripe_rp = stripe_plates[vstripe_idx]
		vstripe = pattern[vstripe_idx]
		_, idx_largest_l_piece_in_vstripe = findmax(bp.l[vstripe])
		# The inner order is not relevant, except the piece with the largest
		# length should be specially treated after. For the same reason we
		# treat the stripe of largest width specially above.
		for (vstripe_idx, piece_idx) in pairs(vstripe)
			vstripe_idx == idx_largest_l_piece_in_vstripe && continue
			plate_idx, stripe_rp = _make_cut!(
				cuts_done, stripe_rp, bp, hcuts_by_pp, bp.l[piece_idx], LENGTH
			)
			_sell_piece!(
				pieces_sold, cuts_done, plate_idx, piece_idx, sllw, bp,
				hcuts_by_pp, vcuts_by_pp, faithful2furini
			)
		end
		# The largest length piece of the largest width stripe can now just be
		# sold. It is cut after the loop because then there is guarantee that all
		# `_make_cut!` calls inside the loop happen in the first half of the
		# current stripe_rp and, therefore, are safe to make. For this last
		# piece a `_make_cut!` call is not necessary and `_sell_piece!` will
		# take care of the edge cases.
		piece_idx = vstripe[idx_largest_l_piece_in_vstripe]
		_sell_piece!(
			pieces_sold, cuts_done, stripe_rp, piece_idx, sllw, bp,
			hcuts_by_pp, vcuts_by_pp, faithful2furini
		)
	end

	return cuts_done, pieces_sold
end
export cuts_and_extractions_from_2_staged_solution

# TODO: document this
function raw_warm_start!(model, nzpe_idxs, nzpe_vals, nzcm_idxs, nzcm_vals)
	@assert length(nzpe_idxs) == length(nzpe_vals)
	@assert length(nzcm_idxs) == length(nzcm_vals)
	pe = model[:picuts]
	cm = model[:cuts_made]
	set_start_value.(pe[nzpe_idxs], nzpe_vals)
	set_start_value.(cm[nzcm_idxs], nzcm_vals)
end
export raw_warm_start!

