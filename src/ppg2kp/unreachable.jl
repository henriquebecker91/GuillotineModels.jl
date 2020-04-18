# highest_plate_idx should be greater-than-or-equal-to the largest
# value found in any field of any element of `cuts`.
function _reachable_plate_types(
	cuts :: Vector{NTuple{3, P}}, highest_plate_idx :: P
) where {P}
	children = [P[] for _ in 1:highest_plate_idx]
	for (pp, fc, sc) in cuts
		push!(children[pp], fc)
		!iszero(sc) && push!(children[pp], sc)
	end
	reached = falses(highest_plate_idx)
	num_reached = zero(P)
	visit_list = P[one(P)]
	next_of_list = one(P)
	while next_of_list <= length(visit_list)
		pp = visit_list[next_of_list]
		if !reached[pp]
			num_reached += one(P)
			reached[pp] = true
			append!(visit_list, [c for c in children[pp] if !reached[c]])
		end
		next_of_list += one(P)
	end
	return num_reached, reached
end

# Works for both cuts and extractions, both have the parent plate as their
# first field.
function _reachable_carves(carves, reachable_plates)
	# We start from the premise the extractions are reachable because this is
	# more probable than the opposite.
	qt_rc = length(carves) # TODO: get the correct number type for this?
	rc = trues(qt_rc)
	for (carve_idx, carve) in pairs(carves)
		pp = first(carve) # gets the parent plate
		if !reachable_plates[pp]
			qt_rc -= one(qt_rc)
			rc[carve_idx] = false
		end
	end
	return qt_rc, rc
end

function _reachable(
	np :: Vector{Tuple{P, D}},
	cuts :: Vector{NTuple{3, P}},
	highest_plate_idx :: P
) where {D, P}
	qt_rp, rp = _reachable_plate_types(cuts, highest_plate_idx)
	qt_rc, rc = _reachable_carves(cuts, rp)
	qt_re, re = _reachable_carves(np, rp)
	return qt_re, re, qt_rc, rc, qt_rp, rp
end

function _print_unreachable_plates(num_reached, reachable, pli_lwb)
	total_num_plates = length(pli_lwb)
	@assert num_reached <= total_num_plates
	@assert length(reachable) == total_num_plates
	qt_unreachable_plate_types = total_num_plates - num_reached
	@show qt_unreachable_plate_types
	println("START_UNREACHABLE_PLATES")
	if qt_unreachable_plate_types > 0
		for (plate_idx, isreachable) in zip(1:total_num_plates, reachable)
			if !isreachable
				l, w, b = pli_lwb[plate_idx]
				println("i $plate_idx l $l w $w b $b")
			end
		end
	end
	println("END_UNREACHABLE_PLATES")
	return
end

function _inner_purge_unreachable!(
	bp :: ByproductPPG2KP{D, S, P}, model, qt_re, re, qt_rc, rc, qt_rp, rp
) where {D, S, P}
	@assert length(bp.np) == length(re)
	@assert length(bp.np) == length(model[:picuts])
	@assert length(bp.cuts) == length(rc)
	@assert length(bp.cuts) == length(model[:cuts_made])
	@assert length(bp.pli_lwb) == length(rp)
	@assert length(bp.pli_lwb) == length(model[:plate_cons])

	# The variables can be deleted first without problem, they make reference
	# to the plates, but the plates have no reference back to the variables.
	qt_re < length(bp.np) && (bp = _delete_vars!(bp, model, re, :picuts))
	qt_rc < length(bp.cuts) && (bp = _delete_vars!(bp, model, rc, :cuts_made))
	# All the code below is just to take care of the case in which some
	# plate types (i.e., constraints) are not reachable (this because all
	# bp.np and bp.cuts elements may need update in this case).
	qt_rp == length(bp.pli_lwb) && return bp

	# NOTE: the code below leaves the positions corresponding to the old index of
	# deleted plates undefined. The assumption is that they will not be accessed,
	# i.e., that no cut or extraction that is kept by this method will ever refer
	# to a removed plate (i.e., that the unreachable plates are, in fact, not
	# reachable from reachable cuts or extractions).
	# `new_plate_idxs[y] == x` where `y` is the old index and `x` the new one.
	new_plate_idxs = Vector{P}(undef, length(bp.pli_lwb))
	next_new_idx = one(P)
	for (old_idx, isreachable) in zip(keys(new_plate_idxs), rp)
		if isreachable
			new_plate_idxs[old_idx] = next_new_idx
			next_new_idx += one(P)
		end
	end

	# Now, for the kept cuts and extractions, we update them to refer to the new
	# id/position of the plates (not the old ones).
	for (cut_idx, (pp, fc, sc)) in pairs(bp.cuts)
		@assert rp[pp]
		@assert rp[fc]
		@assert iszero(sc) || rp[sc]
		new_pp, new_fc = new_plate_idxs[pp], new_plate_idxs[fc]
		new_sc = iszero(sc) ? sc : new_plate_idxs[sc]
		if pp != new_pp || fc != new_fc || sc != new_sc
			bp.cuts[cut_idx] = (new_pp, new_fc, new_sc)
		end
	end
	for (extraction_idx, (pp, ep)) in pairs(bp.np)
		@assert rp[pp]
		new_pp = new_plate_idxs[pp]
		pp != new_pp && (bp.np[extraction_idx] = (new_pp, ep))
	end

	# Then we remove the constraints and the new indexes of cuts and
	# extractions are now correct.
	kept_cons, blot_cons = _partition_by_bits(rp, model[:plate_cons])
	JuMP.delete(model, blot_cons)
	model[:plate_cons] = kept_cons
	deleteat!(bp.pli_lwb, .!rp)

	return bp
end

@timeit TIMER function _purge_unreachable!(bp, model, debug)
	qt_re, re, qt_rc, rc, qt_rp, rp = _reachable(
		bp.np, bp.cuts, lastindex(bp.pli_lwb)
	)
	# TODO: consider printing unreachable vars too?
	debug && _print_unreachable_plates(qt_rp, rp, bp.pli_lwb)
	@assert qt_re <= length(bp.np) && qt_re >= zero(qt_re)
	@assert qt_rc <= length(bp.cuts) && qt_rc >= zero(qt_rc)
	@assert qt_rp <= length(bp.pli_lwb) && qt_rp >= zero(qt_rp)
	if (qt_re + qt_rc + qt_rp) > 0
		if debug
			println("qt_reachable_extractions = $qt_re")
			println("qt_reachable_cuts = $qt_rc")
			println("qt_reachable_plate_types = $qt_rp")
		end
		bp = _inner_purge_unreachable!(bp, model, qt_re, re, qt_rc, rc, qt_rp, rp)
	end
	return bp
end
