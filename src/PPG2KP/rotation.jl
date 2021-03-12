# When working with rotation we have three distinct sets of indexes.
#
# 1. Original piece indexes: they come from the instance and are
#    rotation-unaware (the instance may be used to solve the fixed or
#    rotation-allowed variants).
# 2. Rotation-aware indexes: the indexes of the dummy pieces used to
#    simulate rotation. One original piece may: (i) create just one dummy
#    piece because it is a square or the rotated piece does not fit the
#    large object; (ii) have two dummies but the two dummies are also
#    associated and share demand with another original piece that is
#    the exact rotation of the first original piece; (iii) have two
#    dummies but they are not shared with any other original piece.
#    Consequently the dummy pieces are at least as many as the original
#    pieces and at most the double.
# 3. Shared demand indexes: the indexes of the demand vector in a
#    problem with rotation are always fewer or the same as the number of
#    original (or dummy) pieces. They are smaller if original pieces have
#    their rotation among other original pieces. They may be the same size
#    if all pieces are squares (for example). The smallest this set can go
#    is half the amount of original pieces (in the case every original
#    piece is a rotation of other).
#
#	In the structure below, `rai2sdi` maps (2) to (3), and `sdi2rai` does the
#	inverse (this is, maps (3) to (2)). `shared_demand` is the `d` of a
#	problem with rotation allowed. `sdi2opi` maps (3) to (1), therefore
#	allowing to map (2) to (1) by using `rai2sdi` and then `sdi2opi`.
#	There is no `opi2sdi` just because the code until now has no necessity of
#	this mapping (and it can be deduced from the other if really necessary).
struct RotationAwareData{D, S}
	"Maps rotation-aware indexes to the respective `shared_demand` index."
	rai2sdi :: Vector{D} # rotation-aware piece index to shared demand index
	"Maps `shared_demand` indexes to the respective rotation-aware indexes."
	sdi2rai :: Vector{Union{D, Tuple{D, D}}}
	"Maps `shared_demand` indexes to the respective original piece indexes."
	sdi2opi :: Vector{Union{D, Tuple{D, D}}}
	"The demand of pairs/single rotation-aware pieces."
	shared_demand :: Vector{D}
	# Not used until now, but saved to avoid unnecessary rework.
	"The lengths and widths of the rotation-aware pieces (dummies)."
	ra_sllw :: SortedLinkedLW{D, S}
	# In `get_cut_pattern`, a CutPattern using the dummy pieces is converted to a
	# CutPattern using the original pieces. In the case of perfect rotations
	# within the original pieces, we need to know how much each original piece
	# contributed to the shared demand otherwise the procedure cannot do the
	# conversion in a way that guarantee a valid demand for the original pieces.
	"The demand of the original pieces. Needed for `get_cut_pattern`."
	original_demand :: Vector{D}
end

function _find_perfect_rotation(
	pil :: S, piw :: S, pip :: P, sllw :: SortedLinkedLW{D, S}, p
) :: Int where {D, S, P}
	@assert pil != piw # This code assumes pil x piw is different from piw x pil
	# Search for pieces with an width equal to the piece length and length
	# equal to the piece width (in O(log n) as we are using sorted vectors).
	# (sorted) width/length index range
	swir = searchsorted(sllw.sw, pil)
	slir = searchsorted(sllw.sl, piw)
	# If no piece as pil as width or piw as length then there is no rotation.
	(isempty(swir) || isempty(slir)) && return 0
	# Get the found pieces to the the same numbering (i.e., the instance ordering)
	wir = sllw.swi2pii[swir]
	lir = sllw.sli2pii[slir]
	# Get only the piece indexes that have BOTH pil as width AND piw as length.
	r = intersect(wir, lir)
	# If this set is empty then there is not rotation.
	isempty(r) && return 0
	# Get the profit of the pieces that are rotations.
	pr = p[r]
	# Get the indexes of the right profits in the rotation range.
	cpri = findall(isequal(pip), pr) # correct profit range indexes
	# If no rotation has the same profit, then there are no "perfect rotation"s
	isempty(cpri) && return 0
	# If more than a piece is a perfect rotation then the instance is invalid:
	# there are two pieces with the same dimensions and profits (they should
	# be the same piece with a higher demand).
	if length(cpri) > 1
		idxs = r[cpri]
		@assert allsame(sllw.l[idxs])
		@assert allsame(sllw.w[idxs])
		@assert first(sllw.l[idxs]) == piw
		@assert first(sllw.w[idxs]) == pil
		error(
			"The rotation procedure detected two pieces in the instance that have the same dimensions ($(piw)x$(pil)) and the same profit ($pip). This is considered invalid, both pieces should be aggregated as a single piece with the summed demand. This was only detected because there is a third piece that is a rotation of the other two and this broke some assumptions of the rotation procedure."
		)
	end

	return r[only(cpri)]
end

# If we generalize this to multiple large object (i.e., large objects)
# with distinct dimensions we need to deal with the "rotated small object does
# not fit into the large object" in a better way or just ignore it (i.e.,
# create the dummy rotated object, even if it will not be used).
function build_RAD(
	L :: S, W :: S, sllw :: SortedLinkedLW{D, S}, d :: AbstractVector{D},
	p :: AbstractVector{P}
) :: RotationAwareData{D, S} where {D, S, P}
	@assert length(d) == length(p)
	@assert length(sllw.l) == length(d)

	ral, raw = S[], S[]
	shared_demand = D[]

	rai2sdi = D[]
	sdi2rai = Union{D, Tuple{D, D}}[]
	sdi2opi = Union{D, Tuple{D, D}}[]

	processed = falses(length(d))

	for opi in 1:length(d)
		processed[opi] && continue

		# First, let us add the unrotated dummy piece.
		curr_l = sllw.l[opi]
		curr_w = sllw.w[opi]
		push!(ral, curr_l)
		push!(raw, curr_w)
		# We start assuming the demand of the original piece.
		push!(shared_demand, d[opi])
		push!(rai2sdi, lastindex(shared_demand)) # this cannot change later

		# If the piece is a square or the rotation do not fit, then there
		# is no rotated dummy piece.
		if curr_l == curr_w || curr_l > W || curr_w > L
			push!(sdi2rai, lastindex(rai2sdi))
			push!(sdi2opi, opi)
			processed[opi] = true
			continue
		end

		pr_opi = _find_perfect_rotation(last(ral), last(raw), p[opi], sllw, p)
		if !iszero(pr_opi) # a perfect rotation was found
			# Create the rotated dummy piece. Share its demand.
			push!(ral, curr_w)
			push!(raw, curr_l)
			shared_demand[end] += d[pr_opi]
			push!(rai2sdi, lastindex(shared_demand)) # this one is for the rotation

			# map the demand index to the two dummies, an to the two original
			curr_rai = lastindex(rai2sdi)
			push!(sdi2rai, (curr_rai - 1, curr_rai))
			push!(sdi2opi, (opi, pr_opi))

			# mark both original counterparts as processed
			processed[opi] = true
			processed[pr_opi] = true
			continue
		end

		# If the piece has nothing that hinders rotation, nor another original
		# that is a perfect rotation, then it just spawns a second rotated dummy.
		push!(ral, curr_w)
		push!(raw, curr_l)
		# there is no change to the shared demand here
		push!(rai2sdi, lastindex(shared_demand)) # this one is for the rotation

		# map the demand index to the two dummies, and to the single original
		curr_rai = lastindex(rai2sdi)
		push!(sdi2rai, (curr_rai - 1, curr_rai))
		push!(sdi2opi, opi)

		processed[opi] = true
	end

	return RotationAwareData{D, S}(
		rai2sdi, sdi2rai, sdi2opi, shared_demand, SortedLinkedLW(D, ral, raw),
		deepcopy(d)
	)
end
