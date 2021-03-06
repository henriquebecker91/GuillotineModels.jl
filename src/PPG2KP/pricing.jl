import .Heuristic.fast_iterated_greedy
#import Serialization
#import Gurobi

# Gives the indexes of all elements of bp.cuts that refer to a restricted cut.
# NOTE: no extraction in bp.np can be seen as an unrestricted cut.
function _all_restricted_cuts_idxs(
	bp :: ByproductPPG2KP{D, S, P}
) :: Vector{Int} where {D, S, P}
	rc_idxs = Vector{Int}()
	sizehint!(rc_idxs, length(bp.cuts))
	usl = unique!(sort(bp.l))
	usw = unique!(sort(bp.w))
	n = length(usl)
	for (idx, cuts) in pairs(bp.cuts)
		fc = cuts[FIRST_CHILD]
		if idx >= bp.first_vertical_cut_idx
			# If the cut is vertical it reduces the width.
			if !isempty(searchsorted(usw, bp.pli_lwb[fc][WIDTH]))
				push!(rc_idxs, idx)
			end
		else
			# If the cut is horizontal it reduces the length.
			if !isempty(searchsorted(usl, bp.pli_lwb[fc][LENGTH]))
				push!(rc_idxs, idx)
			end
		end
	end
	return rc_idxs
end

# Internal use.
# Given two iterators, `u` and `s` (both following a common ordering), return
# an iterator over all the elements in `u` but not in `s`. `u` MUST BE
# a superset of `s` or this will not work as intended.
function _setdiff_sorted(u, s)
	isempty(s) && return copy(u)
	t = empty(u)
	next = iterate(s)
	for eu in u
		if next === nothing
			push!(t, eu)
			continue
		end
		(es, ss) = next
		if eu == es
			next = iterate(s, ss)
		else
			push!(t, eu)
		end
	end
	t
end

struct BPLinkage{P}
	exts_part2full :: Vector{P}
	exts_full2part :: Vector{P}
	cuts_part2full :: Vector{P}
	cuts_full2part :: Vector{P}
	plis_part2full :: Vector{P}
	plis_full2part :: Vector{P}
end

#= Use to assert correctness. Uncomment if something seems wrong.
function _check_part2full(part2full)
	@assert allunique(part2full)
	return
end

function _check_full2part(full2part)
	@assert allunique(filter(!iszero, full2part))
	return
end

function _check_linked_idxs(part2full, full2part)
	_check_part2full(part2full)
	_check_full2part(full2part)
	for (part_idx, full_idx) in enumerate(part2full)
		@assert full_idx <= length(full2part)
		@assert full2part[full_idx] == part_idx
	end
	return
end

function _check_linkage(l)
	_check_linked_idxs(l.cuts_part2full, l.cuts_full2part)
	_check_linked_idxs(l.exts_part2full, l.exts_full2part)
	_check_linked_idxs(l.plis_part2full, l.plis_full2part)
	return
end
=#

# Given the type for the index and a bitarray-like, return a duple of linked
# indexes: the first has length equal to the number of true elements in bits
# and the i-th position has the index of the i-th true value in bits; the
# second has length equal to bits, value zero in each position that is false
# in bits, and each position that is true in bits has the number of trues in
# bits up to that position.
function bits2lidxs(::Type{I}, bits) :: NTuple{2, Vector{I}} where {I}
	part2full = I[]
	full2part = Vector{I}(undef, length(bits))
	qt_trues = zero(I)
	for (i, v) in pairs(bits)
		if v
			qt_trues += one(I)
			full2part[i] = qt_trues
			push!(part2full, convert(I, i))
		else
			full2part[i] = zero(I)
		end
	end
	#_check_linked_idxs(part2full, full2part)
	return part2full, full2part
end

@timeit TIMER function _create_restricted_byproduct(
	bp :: ByproductPPG2KP{D, S, P}; copy_unchanged = true,
	given_rc_idxs = nothing, only_reachable = true
) :: Tuple{
	ByproductPPG2KP{D, S, P}, BPLinkage{P}
} where {D, S, P}
	# The cuts/plates/extractions of a restricted model are obtained by a
	# two-step process: (1) get only the cuts in which the cut length/width
	# match the corresponding dimension of a piece; (2) detect which plates,
	# extractions, and restricted cuts are reachable from those cuts.
	# Note that the cuts obtained in (1) can occur in plates that cannot be
	# obtained from the cuts in (1), consequently, the number of restricted cuts
	# often drastically reduces in (2) (and the values become similar to the ones
	# obtained by Furini).
	if isnothing(given_rc_idxs)
		rc_idxs = _all_restricted_cuts_idxs(bp)
	elseif copy_unchanged
		# copyied here because it is returned as a field of BPLinkage
		rc_idxs = copy(given_rc_idxs)
	else
		rc_idxs = given_rc_idxs
	end
	fvc = if bp.first_vertical_cut_idx > length(bp.cuts)
		nothing
	else
		bp.cuts[bp.first_vertical_cut_idx]
	end

	rc_fvci = searchsortedfirst(rc_idxs, bp.first_vertical_cut_idx)
	rc = bp.cuts[rc_idxs]

	if !only_reachable
		kept_fields = (bp.np, bp.pli_lwb, bp.d, bp.l, bp.w, bp.L, bp.W)
		if copy_unchanged
			kept_fields = copy.(kept_fields)
		end
		new_bp = ByproductPPG2KP{D, S, P}(
			rc, rc_fvci, kept_fields..., bp.mirror_plates
		)
		cuts_full2part = zeros(P, length(bp.cuts))
		cuts_full2part[rc_idxs] .= 1:length(rc_idxs)
		exts = collect(1:length(bp.np))
		plis = collect(1:length(bp.pli_lwb))
		bpl = BPLinkage(
			exts, copy(exts), rc_idxs, cuts_full2part, plis, copy(plis)
		)
		return new_bp, bpl
	end

	qt_re, re_bits, qt_rc, rc_bits, qt_rp, rp_bits = _reachable(
		bp.np, rc, lastindex(bp.pli_lwb)
	)
	# The rc_bits cannot be directly used to create the BPLinkage,
	# this because they refer to `rc` not to `bp.cuts` (i.e., a subset
	# not the whole). Some indexing magic is necessary to convert the
	# bitarray of the subset to a bitarray of the whole.
	rc_full_bits = falses(length(bp.cuts))
	rc_full_bits[rc_idxs] .= rc_bits
	# `re_bits` refer to `bp.np` and `rp_bits` refer to
	# `1:length(lastindex(bp.pli_lwb))` so they are fine.
	exts_part2full, exts_full2part = bits2lidxs(P, re_bits)
	cuts_part2full, cuts_full2part = bits2lidxs(P, rc_full_bits)
	plis_part2full, plis_full2part = bits2lidxs(P, rp_bits)
	bpl = BPLinkage(
		exts_part2full, exts_full2part, cuts_part2full, cuts_full2part,
		plis_part2full, plis_full2part
	)
	#_check_linkage(bpl)
	re = bp.np[exts_part2full]
	rp = bp.pli_lwb[plis_part2full]
	deleteat!(rc, .!rc_bits) # rc is already a copy of a subset
	rc_fvci = rc_fvci > length(rc_idxs) ? length(rc) + 1 : sum(rc_bits[1:rc_fvci])
	# If this assert fails, it means our index tricks using searchsortedfirst
	# and sum have failed us.
	@assert rc_fvci > length(rc) || rc[rc_fvci] == fvc
	# After the assert, we can update the values of the cuts and extractions to
	# refer to the new plates.
	# TODO: put these loops inside a internal method of their own, I have
	# the feeling I have already written them many times.
	for (i, (pp, fc, sc)) in pairs(rc)
		new_pp, new_fc = plis_full2part[pp], plis_full2part[fc]
		new_sc = iszero(sc) ? sc : plis_full2part[sc]
		#@assert new_pp > 0 && new_pp <= length(rp)
		#@assert new_fc > 0 && new_fc <= length(rp)
		#@assert new_sc >= 0 && new_fc <= length(rp)
		rc[i] = (new_pp, new_fc, new_sc)
	end
	for (i, (plate, piece)) in pairs(re)
		#@assert !iszero(oldpli2newpli[plate])
		re[i] = (plis_full2part[plate], piece)
		#@assert re[i][1] > 0 && re[i][1] <= length(rp)
	end

	kept_fields = (bp.d, bp.l, bp.w, bp.L, bp.W)
	copy_unchanged && (kept_fields = copy.(kept_fields))
	r_bp = ByproductPPG2KP{D, S, P}(
		rc, rc_fvci, re, rp, kept_fields..., bp.mirror_plates
	)
	return r_bp, bpl
end

# Arguments:
# * `cut`: a triple of `pp`, `fc`, and `sc`.
#   + `pp`: parent plate id/index.
#   + `fc`: first child id/index.
#   + `sc`: second child id/index (zero means there is no second child).
# * `constraints`: the constraints related to the plates; if indexed by
#   the plate id, it gives the corresponding constraint.
function _reduced_profit(
	cut :: Tuple{P, P, P}, constraints :: Vector{T}
) :: Float64 where {P, T}
	pp, fc, sc = cut
	pp_dual = dual(constraints[pp])
	fc_dual = dual(constraints[fc])
	sc_dual :: Float64 = iszero(sc) ? 0.0 : dual(constraints[sc])
	# The reduced profit computation shown in the paper is done here.
	rp = -fc_dual + -sc_dual - -pp_dual
	#=if rp < zero(rp)
		@show cut
		@show pp_dual
		@show fc_dual
		@show sc_dual
		@show rp
	end=#
	return rp
end

function _reduced_profit(
	cut :: Tuple{P, P, P}, con_duals :: Vector{Float64}
) :: Float64 where {P}
	pp, fc, sc = cut
	pp_dual = con_duals[pp]
	fc_dual = con_duals[fc]
	sc_dual :: Float64 = iszero(sc) ? 0.0 : con_duals[sc]
	# The reduced profit computation shown in the paper is done here.
	rp = -fc_dual + -sc_dual - -pp_dual
	return rp
end

@timeit TIMER function _restricted_final_pricing!(
	model, p :: Vector{P},
	bp :: ByproductPPG2KP{D, S, P}, lidxs :: BPLinkage{P},
	seed, verbose :: Bool, mip_start :: Bool,
	bm :: BaseModel, start :: Float64 = time(),
	limit :: Float64 = float(60*60*24*365)
) where {D, S, P}
	bp = deepcopy(bp) # avoids tainting this parameter
	n = num_variables(model)
	pe = model[:picuts] # Piece Extractions
	cm = model[:cuts_made] # Cuts Made
	# This assert is here because, if this becomes false some day, the code
	# will need to be reworked. Now it only relax/fix variables in those sets
	# so it will need to be sensibly extended to new sets of variables.
	@assert n == length(cm) + length(pe)
	# Relax all variables.
	@timeit TIMER "relax_cm" cm_svcs = relax!(cm)
	@timeit TIMER "relax_pe" pe_svcs = relax!(pe)
	@assert length(cm) == length(bp.cuts)

	# The 'best known value'/'primal bound' of the heuristic is needed to do
	# the pricing. If we can MIP-start the model, we can just call
	# `mip_start_by_heuristic!` and get the `bkv` returned, otherwise we need
	# to call just the heuristic to get the bkv (and throw away the rest).
	heuristic_lb_time = @elapsed if mip_start
		(bkv, _, _), heuristic_raw_ws = mip_start_by_heuristic!(
			model, bp, p, seed, bm
		)
	else
		bkv :: P, _, _ = fast_iterated_greedy(
			bp.d, p, bp.l, bp.w, bp.L, bp.W, Xoroshiro128Plus(seed)
		)
	end
	restricted_LB = LB = convert(Float64, bkv)
	if verbose
		@show restricted_LB
		println("heuristic_lb = $(restricted_LB)")
		@show heuristic_lb_time
	end

	# Check if the heuristic did not blow the time limit, but only after
	# the information output that happens after it.
	throw_if_timeout_now(start, limit)

	# Solve the relaxed restricted model.
	verbose && println("MARK_FURINI_PRICING_RESTRICTED_LP_SOLVE")
	@timeit TIMER "lp_solve" optimize_within_time_limit!(
		model, limit - (time() - start)
	)
	if verbose
		restricted_lp_stop_reason = termination_status(model)
		@show restricted_lp_stop_reason
		restricted_lp_stop_code = Int(restricted_lp_stop_reason)
		@show restricted_lp_stop_code
	end
	if primal_status(model) == MOI.FEASIBLE_POINT
		restricted_UB = LP = objective_value(model)
		verbose && @show restricted_UB
	end

	throw_if_timeout_now(start, limit)

	# Check if everything seems ok with the values obtained.
	if termination_status(model) != MOI.OPTIMAL
		error(
			"For some reason, solving the relaxed restricted model did not" *
			" terminate with optimal status (the status was" *
			" $(termination_status(model))). The code is not prepared to deal" *
			" with this possibility and will abort."
		)
	end
	# Note: the iterated_greedy solve the shelf version of the problem,
	# that is restricted, so it cannot return a better known value than the
	# restricted "upper bound"/"linear(ized) problem"/"integer relaxation".
	@assert restricted_LB <= restricted_UB + eps(restricted_UB)
	# This is a hack. If restricted_LB (that is guaranteed to be an integer
	# considering our assumption that piece profits are integer) is less than
	# one unity smaller than restricted_UB, then the optimal for the
	# restricted problem was already found, and we do not need to price the
	# restricted model and solve it.
	if restricted_UB - restricted_LB < 1.0
		# Keep the MIP-start in full indexes, so we can just convert it back to
		# the final indexes, instead of changing each time the part indexes change.
		heuristic_raw_ws[1] .= lidxs.exts_part2full[heuristic_raw_ws[1]]
		heuristic_raw_ws[3] .= lidxs.cuts_part2full[heuristic_raw_ws[3]]
		# If LB and UB match, the restricted priced model is an empty model.
		return bkv, heuristic_raw_ws, true #= early_return =#, bp
	end

	plate_duals = dual.(model[:plate_cons])
	kept, bp = _final_pricing!(model, plate_duals, bp, LB, LP)

	# Update the BPLinkage (some cuts were deleted). Needed to convert the
	# MIP-start part indexes to full indexes before returning.
	_delete_from_part!(lidxs.cuts_part2full, lidxs.cuts_full2part, .!kept)

	# Already prints this info before we risk timeout solving it.
	if verbose
		qt_cmvars_deleted_by_heuristic_pricing = length(kept) - sum(kept)
		@show qt_cmvars_deleted_by_heuristic_pricing
		println("qt_cmvars_priced_restricted = $(length(bp.cuts))")
		println("qt_pevars_priced_restricted = $(length(bp.np))")
		println("qt_plates_priced_restricted = $(length(bp.pli_lwb))")
	end

	# Make the model a MIP again.
	restore!.(model[:picuts], pe_svcs)
	restore!.(model[:cuts_made], cm_svcs[kept])

	# restricted MIP solved
	verbose && println("MARK_FURINI_PRICING_RESTRICTED_MIP_SOLVE")
	@timeit TIMER "mip_solve" optimize_within_time_limit!(
		model, limit - (time() - start)
	)
	# If we MIP start the restricted model with a feasible solution, it should
	# be impossible to get a different status here.
	@assert !mip_start || primal_status(model) == MOI.FEASIBLE_POINT
	# The assert above is valid, but this does not mean we actually deal
	# with the possibility of not having a primal solution below.
	model_obj = round(P, objective_value(model), RoundNearest)
	if verbose
		restricted_stop_reason = termination_status(model)
		restricted_stop_code = Int(restricted_stop_reason)
		println("restricted_obj_value = $model_obj")
		println("restricted_obj_bound = $(objective_bound(model))")
		@show restricted_stop_reason
		@show restricted_stop_code
	end

	throw_if_timeout_now(start, limit)

	model_obj > bkv && (bkv = model_obj)
	# Save the restricted model optimal solution so the final model can be
	# warm-start or, if iterative pricing proves this solution is optimal,
	# return this solution.
	restricted_sol = save_mip_start(model)
	# Convert it to full indexes before returning.
	restricted_sol[1] .= lidxs.exts_part2full[restricted_sol[1]]
	restricted_sol[3] .= lidxs.cuts_part2full[restricted_sol[3]]

	return bkv, restricted_sol, false #= early_return =#, bp
end

# Arguments:
# * `pool_idxs_to_add`: object that may or not be empty, but will have its
#   old content thrown away, be resized, and in the end will have all pool
#   indexes corresponding to cuts/vars that should be added to the model.
# * `pool`: pool of all cuts that are outside the model.
# * `plate_cons`: vector of constraints representing the plates in the model.
# * `threshold`: reduced profit threshold used in the process.
# * `n_max`: maximum size of pool_idxs_to_add (i.e., maximum amount of
#   variable indexes that we are interested in adding in a single step).
# The most relevant byproduct of this method are the changes to
# pool_idxs_to_add (i.e., which vars need to be unfixed).
# Return: the number of positive reduced profit variables found,
# and a boolean indicating if some variable above the threshold was found.
@timeit TIMER function _recompute_idxs_to_add!(
	pool_idxs_to_add, pool, plate_duals, threshold, n_max :: P,
	sort_by_rp :: Bool, rps_buffer = Vector{Float64}(undef, length(pool))
) :: Tuple{P, Bool} where {P}
	pool_size = length(pool)
	if sort_by_rp && length(rps_buffer) < pool_size
		resize!(rps_buffer, pool_size)
	end
	# Variables used in the loop.
	found_above_threshold = false
	num_positive_rp_vars = zero(P)
	qt_selected_vars = zero(P)
	for (pool_idx, cut) in pairs(pool)
		rp = _reduced_profit(cut, plate_duals)
		sort_by_rp && (rps_buffer[pool_idx] = rp)
		# non-positive reduced profit is irrelevant to us
		rp <= 0.0 && continue
		num_positive_rp_vars += one(P)
		if rp > threshold
			if !found_above_threshold
				qt_selected_vars = zero(P)
				found_above_threshold = true
			end
			qt_selected_vars += one(P)
			if qt_selected_vars > length(pool_idxs_to_add)
				push!(pool_idxs_to_add, pool_idx)
			else
				pool_idxs_to_add[qt_selected_vars] = pool_idx
			end
			# If n_max variables above the threshold exist, and we do not
			# sort by reduced profit, only them are used.
			qt_selected_vars >= n_max && !sort_by_rp && break
		elseif !found_above_threshold
			# Unfortunately we cannot break here if we find n_max variables because
			# we can find one above the threshold later yet. At least, we can stop
			# pushing new variables to the pool.
			qt_selected_vars >= n_max && !sort_by_rp && continue
			qt_selected_vars += one(P)
			if qt_selected_vars > length(pool_idxs_to_add)
				push!(pool_idxs_to_add, pool_idx)
			else
				pool_idxs_to_add[qt_selected_vars] = pool_idx
			end
		end
	end
	#vector_summary(positive_rps)

	return_size = min(qt_selected_vars, n_max)
	sort_by_rp && partialsort!(
		(@view pool_idxs_to_add[1:qt_selected_vars]), 1:return_size;
		rev = true, by = i -> rps_buffer[i]
	)
	# Only n_max variables are added/unfixed in a single iteration.
	resize!(pool_idxs_to_add, return_size)
	# deleteat! needs the indexes to be sorted
	sort_by_rp && sort!(pool_idxs_to_add)

	return num_positive_rp_vars, found_above_threshold
end

@timeit TIMER function _safe_recompute_idxs_to_add!(
	pool, plate_duals, threshold, n_max :: P, sort_by_rp :: Bool
) :: Tuple{Vector{P}, P, Bool} where {P}
	pool_size = length(pool)
	sort_by_rp && (rps_buffer = Vector{Float64}(undef, pool_size))
	pool_idxs_to_add = P[]
	# Variables used in the loop.
	found_above_threshold = false
	num_positive_rp_vars = zero(P)
	qt_selected_vars = zero(P)
	for (pool_idx, cut) in pairs(pool)
		rp = _reduced_profit(cut, plate_duals)
		sort_by_rp && (rps_buffer[pool_idx] = rp)
		# non-positive reduced profit is irrelevant to us
		rp <= 0.0 && continue
		num_positive_rp_vars += one(P)
		if rp > threshold
			if !found_above_threshold
				qt_selected_vars = zero(P)
				found_above_threshold = true
			end
			push!(pool_idxs_to_add, pool_idx)
			# If n_max variables above the threshold exist, and we do not
			# sort by reduced profit, only them are used.
			length(pool_idxs_to_add) >= n_max && !sort_by_rp && break
		elseif !found_above_threshold
			# Unfortunately we cannot break here if we find n_max variables because
			# we can find one above the threshold later yet. At least, we can stop
			# pushing new variables to the pool.
			length(pool_idxs_to_add) >= n_max && !sort_by_rp && continue
			push!(pool_idxs_to_add, pool_idx)
		end
	end
	#vector_summary(positive_rps)

	return_size = min(length(pool_idxs_to_add), n_max)
	sort_by_rp && partialsort!(pool_idxs_to_add, 1:return_size;
		rev = true, by = i -> rps_buffer[i]
	)
	# Only n_max variables are added/unfixed in a single iteration.
	resize!(pool_idxs_to_add, return_size)
	# deleteat! needs the indexes to be sorted
	sort_by_rp && sort!(pool_idxs_to_add)

	return pool_idxs_to_add, num_positive_rp_vars, found_above_threshold
end

# The values in `to_delete` refer to the indexes of `part2full`
# (not the values of `part2full` that are indexes of `full2part`).
function _delete_from_part!(part2full, full2part, to_delete)
	# Supports BitArray and Vector{Bool} by creating an intermediary array, maybe
	# optimize to have multiple method definitions in the future.
	idxs_to_delete = keys(part2full)[to_delete]
	isempty(idxs_to_delete) && return part2full
	# Quantity of deleted indexes before the current index.
	qt_to_shift = zero(eltype(full2part))
	# The loop acts in the ranges between the first deleted index and the
	# second deleted index, the second and the thirds, and so on. The
	# undeleted values in such ranges need full2part to be updated so
	# it points to the right positions after the deletion.
	for del_list_idx = 1:(length(idxs_to_delete) - 1)
		deleted_idx = idxs_to_delete[del_list_idx]
		full2part[part2full[deleted_idx]] = zero(eltype(full2part))
		qt_to_shift += one(qt_to_shift)
		range_start = deleted_idx + 1
		range_end = idxs_to_delete[del_list_idx + 1] - 1
		for part_idx in range_start:range_end
			full2part[part2full[part_idx]] -= qt_to_shift
		end
	end
	# The loop above does not take care of the range from the last deleted index
	# to the end of part2full. The code below takes care of it.
	deleted_idx = last(idxs_to_delete)
	full2part[part2full[deleted_idx]] = zero(eltype(full2part))
	qt_to_shift += one(qt_to_shift)
	@assert qt_to_shift == length(idxs_to_delete)
	range_start = deleted_idx + 1
	range_end = lastindex(part2full)
	for part_idx in range_start:range_end
		full2part[part2full[part_idx]] -= qt_to_shift
	end
	# Finally, the values in part2full are really deleted.
	deleteat!(part2full, to_delete)
	#_check_linked_idxs(part2full, full2part)
	return part2full
end

# The values in `to_append` refer to values in `part2full` (i.e., indexes in
# `full2part`). `to_append` should contain no values already found in
# `part2full` before the append.
function _append_to_part!(part2full, full2part, to_append)
	l = length(part2full)
	for (i, v) in enumerate(to_append)
		full2part[v] = l + i
	end
	append!(part2full, to_append)
	#_check_linked_idxs(part2full, full2part)
	return part2full
end

# Given a PPG2KP model that has `cut_var`, and in which all constraints that
# `cut` references are present in `plate_cons`, set the variable coefficient in
# the referenced constraints to the right value for a PPG2KP model.
function _set_cut_coeffs!(cut, cut_var, plate_cons)
	pp, fc, sc = cut
	set_normalized_coefficient(plate_cons[pp], cut_var, 1.0)
	if iszero(sc)
		set_normalized_coefficient(plate_cons[fc], cut_var, -1.0)
	elseif fc == sc
		# This is the case in which the plate is divides by the exact
		# middlepoint and births two copies of the same child plate.
		set_normalized_coefficient(plate_cons[fc], cut_var, -2.0)
	else
		set_normalized_coefficient(plate_cons[fc], cut_var, -1.0)
		set_normalized_coefficient(plate_cons[sc], cut_var, -1.0)
	end
	return
end

function _set_extraction_coeffs!(extraction, var, plate_cons, demand_cons)
	plate, piece = extraction
	set_normalized_coefficient(plate_cons[plate], var, 1.0)
	set_normalized_coefficient(demand_cons[piece], var, 1.0)
end

# Given a PPG2KP model (with a :plate_cons structure associated)
# and a number of new plate_cons to add, the method creates
# `qt_plates_cons_to_add` empty constraints and append them to
# model[:plate_cons] while setting their name correctly.
function _add_plate_cons!(
	model, qt_new_cons
) :: Nothing
	iszero(qt_new_cons) && return
	plate_cons = model[:plate_cons]
	qt_old_cons = length(plate_cons)
	new_cons = @constraint(model, [i=1:qt_new_cons], 0.0 <= 0.0)
	new_cons_idxs = (qt_old_cons + 1):(qt_old_cons + qt_new_cons)
	@. set_name(new_cons, "plate_cons[" * string(new_cons_idxs) * "]")
	append!(plate_cons, new_cons)
	return
end

# NOTE: pricing may not work with hybridization because of the below
# (i.e., it assumes non-zero fc).
# Pattern repeated in many parts of the code.
function _map_cuts(f, cuts)
	return map(cuts) do (pp, fc, sc)
		(f(pp), f(fc), iszero(sc) ? sc : f(sc))
	end
end

function _add_cuts_from_pool!(
	model, full_cut_idxs_to_add :: AbstractVector{P},
	pli2pair :: Vector{Vector{P}}, lidxs :: BPLinkage{P},
	full_bp :: ByproductPPG2KP{D, S, P}, part_bp :: ByproductPPG2KP{D, S, P}
) where {D, S, P}
	# Shorten the names.
	picuts = model[:picuts]
	plate_cons = model[:plate_cons]
	demand_con = model[:plate_cons]

	# Get some basic info from the parameters.
	qt_new_cuts = length(full_cut_idxs_to_add)
	cuts_to_add = full_bp.cuts[full_cut_idxs_to_add]

	# ==================== PLATE CONSTRAINTS ====================

	# Before anything else, we check if these cuts create plates that did not
	# exist before in the part(ial) model.

	# Get the full indexes of the that plates will exist in the model
	# only when these cuts were added.
	full_plis_to_add = P[]
	@unpack plis_part2full, plis_full2part = lidxs
	# (f)ull {(p)arent (p)late, {(f)irst, (s)econd} (c)hild}
	for (fpp, ffc, fsc) in cuts_to_add # has the full indexes
		@assert !iszero(ffc) && !iszero(fpp)
		iszero(plis_full2part[fpp]) && push!(full_plis_to_add, fpp)
		iszero(plis_full2part[ffc]) && push!(full_plis_to_add, ffc)
		iszero(fsc) || iszero(plis_full2part[fsc]) && push!(full_plis_to_add, fsc)
	end
	# The same plates not yet in the model can appear in multiple cuts,
	# so the vector needs to be cleaned up so each constraint is added a
	# single time.
	unique!(sort!(full_plis_to_add))

	# Create the missing plate constraints, update their linked indexes,
	# add the plates to part_bp.pli_lwb.
	_add_plate_cons!(model, length(full_plis_to_add))
	_append_to_part!(plis_part2full, plis_full2part, full_plis_to_add)
	append!(part_bp.pli_lwb, full_bp.pli_lwb[full_plis_to_add])

	# ==================== CUTS ====================

	# Create the new cut variables and add them to the model.
	qt_cuts_before = length(lidxs.cuts_part2full)
	upper_bounds = getindex.(full_bp.pli_lwb[getindex.(cuts_to_add, PARENT)], 3)
	new_cut_vars = @variable(
		model, [i=1:qt_new_cuts], lower_bound = 0.0, upper_bound = upper_bounds[i]
	)
	new_cut_vars_range = (qt_cuts_before + 1):(qt_cuts_before + qt_new_cuts)
	@. set_name(new_cut_vars, "cuts_made[" * string(new_cut_vars_range) * "]")
	cuts_made = model[:cuts_made]
	append!(cuts_made, new_cut_vars)

	# Update the cut index linking between full and part models.
	@unpack cuts_part2full, cuts_full2part = lidxs
	_append_to_part!(cuts_part2full, cuts_full2part, full_cut_idxs_to_add)

	# Create the cuts in the partial model, they are not guaranteed to be the
	# same objects, as they are just triples of plate indexes, and the same
	# plate can have different indexes in the full and part(ial) models.
	cuts_to_add_in_part = _map_cuts(cuts_to_add) do pli
		lidxs.plis_full2part[pli]
	end
	append!(part_bp.cuts, cuts_to_add_in_part)

	# Use the cuts with corrected indexes to set the coefficients of the
	# newly created cut variables in the newly created plate constraints.
	_set_cut_coeffs!.(cuts_to_add_in_part, new_cut_vars, (plate_cons,))

	# ==================== EXTRACTIONS ====================

	# Discover which extractions do not exist yet in the (part)ial model.
	new_full_exts_idxs = collect(Iterators.flatten(pli2pair[full_plis_to_add]))
	new_exts = full_bp.np[new_full_exts_idxs]
	qt_new_exts = length(new_exts)

	# If there is new extraction variables to be added...
	if !iszero(qt_new_exts)
		@unpack exts_part2full, exts_full2part = lidxs

		# Add the extraction vars to the model.
		upper_bounds = map(new_full_exts_idxs) do i
			# The field `1` in a `np` element is the plate index, the `3` in
			# a `pli_lwb` element is the (b)ound. The field `2` in a np element
			# is the piece. TODO: someday make it readable.
			min(full_bp.pli_lwb[full_bp.np[i][1]][3], full_bp.d[full_bp.np[i][2]])
		end
		qt_exts_before = length(exts_part2full)
		new_exts_vars = @variable(
			model, [i=1:qt_new_exts], lower_bound = 0.0,
			upper_bound = upper_bounds[i]
		)
		new_exts_vars_range = (qt_exts_before + 1):(qt_exts_before + qt_new_exts)
		@. set_name(new_exts_vars, "picuts[" * string(new_exts_vars_range) * "]")
		append!(picuts, new_exts_vars)

		# Update the Linkage of the extraction vars. Needed for next step.
		_append_to_part!(exts_part2full, exts_full2part, new_full_exts_idxs)

		# Set the coefficient of the extraction vars in the constraints.
		# Convert the full plate indexes inside the extractions to the indexes in
		# the partial model, use these "partialized" extractions to set the
		# coefficients in the partial model.
		exts_to_add_in_part = map(new_exts) do (plate, piece)
			(lidxs.plis_full2part[plate], piece)
		end
		_set_extraction_coeffs!.(
			exts_to_add_in_part, new_exts_vars, (plate_cons,), (demand_con,)
		)
		# Add the extractions to the part_bp.
		append!(part_bp.np, exts_to_add_in_part)
	end

	return
end

# Given the linked indexes part2full and full2part, sort the values of
# part2full (i.e., stored indexes of full2part) in increasing order,
# update full2part to point to the correct positions, and sort all
# vectors in part_vecs the same way. Consequently, the subset of the full
# vectors represented by part2full and part_vecs now has the same
# order of their full counterparts.
function _sort_linkage!(part2full, full2part, part_vecs)
	sp = sortperm(part2full)
	position_when_sorted = invperm(sp)
	for (part_idx, full_idx) in enumerate(part2full)
		full2part[full_idx] = position_when_sorted[part_idx]
	end
	@assert issorted(filter(!iszero, full2part))
	part2full .= part2full[sp]
	for part_vec in part_vecs
		part_vec .= part_vec[sp]
	end

	#_check_linked_idxs(part2full, full2part)
	return
end

# Create a new byproduct to replace the parameter part_bp, the
# model (in truth, the variable and constraint vectors associated)
# and the linked indexes are updated to refer to this new ByproductPPG2KP.
# The `full_bp` is not changed, and `part_bp` is invalidated (its vectors
# are reused, so it should not be used anymore after this).
function _clean_partial_byproduct!(model, lidxs, full_bp, part_bp)
	# The part_bp.cuts and part_bp.np are special in the context of this
	# method, they have references to the plate indexes, so their value (not
	# only position) needs to be changed (as the plate indexes will change).
	oldidx2newidx = invperm(sortperm(lidxs.plis_part2full))
	part_bp.cuts .= _map_cuts(i -> oldidx2newidx[i], part_bp.cuts)
	part_bp.np .= map(part_bp.np) do (pli, pii)
		(oldidx2newidx[pli], pii)
	end
	_sort_linkage!(
		lidxs.cuts_part2full, lidxs.cuts_full2part,
		(part_bp.cuts, model[:cuts_made])
	)
	_sort_linkage!(
		lidxs.exts_part2full, lidxs.exts_full2part,
		(part_bp.np, model[:picuts])
	)
	_sort_linkage!(
		lidxs.plis_part2full, lidxs.plis_full2part,
		(part_bp.pli_lwb, model[:plate_cons])
	)
	#_check_linkage(lidxs)
	new_fvci = searchsortedfirst(
		lidxs.cuts_part2full, full_bp.first_vertical_cut_idx
	)
	return ByproductPPG2KP( # new_fvci was computed and need to be updated
		part_bp.cuts, new_fvci, part_bp.np, part_bp.pli_lwb, part_bp.d,
		part_bp.l, part_bp.w, part_bp.L, part_bp.W, part_bp.mirror_plates
	)
end

# Structs and functions to store and print stats related to the model.
#=
abstract type AbstractAggregate{T} end;
function Base.getindex(a :: AbstractAggregate{T}, s :: Symbol) :: T where {T}
	return (getfield(a, s) :: T)
end
struct ModelStats{T} <: AbstractAggregate{T}
	cuts :: T
	extractions :: T
	plates :: T
end
=#

const ModelStats{P} = NamedTuple{(:cmvars, :pevars, :plates), NTuple{3, P}}
function ModelStats(
	bp :: ByproductPPG2KP{D, S, P}
) :: ModelStats{P} where {D, S, P}
	return ModelStats{P}(convert.(P, length.((bp.cuts, bp.np, bp.pli_lwb))))
end
mutable struct IterativePricingStats{P}
	starting_qt :: ModelStats{P}
	previous_qt :: ModelStats{P}
	current_iter :: P
	qt_threshold_hits :: P
end
function IterativePricingStats(
	bp :: ByproductPPG2KP{D, S, P}
) :: IterativePricingStats{P} where {D, S, P}
  initial_stats = ModelStats(bp)
	return IterativePricingStats(initial_stats, initial_stats, zero(P), zero(P))
end
function _print_and_update!(
	stats :: IterativePricingStats{P},
	bp :: ByproductPPG2KP{D, S, P},
	threshold_hit :: Bool
) :: Nothing where {D, S, P}
	println("threshold_hit_in_iter_$(stats.current_iter) = $(threshold_hit)")
	stats.qt_threshold_hits += threshold_hit
	curr_stats = ModelStats(bp)
	for name in fieldnames(ModelStats{P})
		println(
			"qt_$(name)_added_by_iter_$(stats.current_iter) = " *
			string(curr_stats[name] - stats.previous_qt[name])
		)
		println(
			"qt_$(name)_after_iter_$(stats.current_iter) = " *
			string(curr_stats[name])
		)
	end
	stats.previous_qt = curr_stats
	stats.current_iter += 1

	return
end

# NOTE: expect the model to already be relaxed, and with only the restricted
# cut variables not fixed to zero (while the rest is fixed to zero).
# NOTE: the description of this method (Explained in 10.1287/ijoc.2016.0710,
# p. 13 (747), last paragraph before section 4.3.) is not entirely clear.
# If n_max is larger than the number of variables above the threshold, then
# should only the variables above the threshold to be added, or should be
# used variables below the threshold to complete the quantity defined by n_max?
# Also, the paper mention the "first X variables", so it is not worth getting
# all variables and ordering them by reduced profit (if they are more than
# n-max)?
@timeit TIMER function _iterative_pricing!(
	model, full_bp :: ByproductPPG2KP{D, S, P},
	part_bp :: ByproductPPG2KP{D, S, P}, lidxs :: BPLinkage{P},
	full_pli2pair :: Vector{Vector{P}}, max_profit :: P,
	alpha :: Float64, beta :: Float64, sort_by_rp :: Bool, debug :: Bool = false,
	start :: Float64 = time(), limit :: Float64 = float(60*60*24*365)
) :: ByproductPPG2KP{D, S, P} where {D, S, P}
	part_bp = deepcopy(part_bp) # Avoid invalidating the original part_bp.
	stats = IterativePricingStats(part_bp)

	flush_all_output()
	# the last solve before this was MIP and has no duals
	debug && println("MARK_FURINI_PRICING_ITERATED_LP_SOLVE_0")
	@timeit TIMER "solve_lp" optimize_within_time_limit!(model, start, limit)
	flush_all_output()

	# The p̄ value in the paper, if there are variables with reduced profit
	# above this threshold then they are added (and none below the threshold),
	# but the n_max limit of variables added in a single iteration is respected.
  threshold = max_profit * beta
	@assert threshold > zero(threshold) # threshold is always larger than zero

	# Gets just the cuts of the full model that are absent in the partial model.
	@assert issorted(lidxs.cuts_part2full)
	unused_full_cuts_idxs = _setdiff_sorted(
		keys(full_bp.cuts), lidxs.cuts_part2full
	)
	unused_full_cuts = full_bp.cuts[unused_full_cuts_idxs]

	plate_cons = model[:plate_cons]
	plate_duals = zeros(Float64, length(full_bp.pli_lwb))
	plate_duals[lidxs.plis_part2full] .= dual.(plate_cons)

	# Do the initial pricing, necessary to compute n_max, and that is always done
	# (i.e., the end condition can only be tested after this first loop).
	pool_idxs_to_add, initial_num_positive_rp_vars, was_above_threshold =
		_safe_recompute_idxs_to_add!(
			unused_full_cuts, plate_duals, threshold, typemax(P), sort_by_rp
		)

	# NOTE: the _recompute_idxs_to_add! below was deemed excessive optimization.
	# The vectors are all allocated here (but used inside
	# `_recompute_idxs_to_add!`) to avoid reallocating them every loop.
	# Note that the `pool_idxs_to_add` values are indexes for the
	# `unused_full_cuts_idxs` vector, they are not the cut indexes themselves.
	#pool_idxs_to_add = P[]
	#rps_buffer = Vector{Float64}(undef, length(unused_full_cuts))
	#initial_num_positive_rp_vars, was_above_threshold = _recompute_idxs_to_add!(
	#	pool_idxs_to_add, unused_full_cuts, plate_duals, threshold, typemax(P),
	#	sort_by_rp, rps_buffer
	#)

	debug && @show initial_num_positive_rp_vars
	# The maximum number of variables unfixed at each iteration.
	n_max = round(P, initial_num_positive_rp_vars * alpha, RoundUp)
	debug && @show n_max
	# As we called `num_positive_rp_vars!` with `typemax(P)` instead `n_max`
	# to be able to compute `n_max` in the first place, we need to resize
	# this first list to `n_max` (in the loop below `_recompute_idxs_to_add!`
	# will do it for us).
	n_max < length(pool_idxs_to_add) && resize!(pool_idxs_to_add, n_max)
	num_positive_rp_vars = initial_num_positive_rp_vars
	# The iterative pricing continue until there are variables to unfix.
	while !isempty(pool_idxs_to_add)
		# Updates the model, part_bp, and lidxs. The part_bp is left in a
		# not-quite-valid-state, vertical and horizontal cuts intercalate and the
		# value of first_vertical_cut_idx is bullshit. But such state is enough for
		# passing it again to _add_cuts_from_pool. It is fixed just before return.
		_add_cuts_from_pool!(
			model, unused_full_cuts_idxs[pool_idxs_to_add], full_pli2pair, lidxs,
			full_bp, part_bp
		)
		if debug
			_print_and_update!(stats, part_bp, was_above_threshold)
			println(
				"qt_positive_rp_cmvars_iter_$(stats.current_iter) = " *
				string(num_positive_rp_vars)
			)
		end

		#_check_linkage(lidxs)

		# Now that the cuts were added to the model from the pool they can be
		# removed from the pool.
		@timeit TIMER "deleteat!" begin
			deleteat!(unused_full_cuts_idxs, pool_idxs_to_add)
			deleteat!(unused_full_cuts, pool_idxs_to_add)
		end
		@assert length(unused_full_cuts_idxs) == length(unused_full_cuts)

		# Solve the LP model again so the constraint duals are updated.
		flush_all_output()
		if debug
			println("MARK_FURINI_PRICING_ITERATED_LP_SOLVE_$(stats.current_iter)")
		end

		@timeit TIMER "solve_lp" optimize_within_time_limit!(model, start, limit)
		flush_all_output()

		# The only information needed from the model solving.
		plate_duals[lidxs.plis_part2full] .= dual.(plate_cons)

		# Excessively optimized version disabled.
		#num_positive_rp_vars, was_above_threshold = _recompute_idxs_to_add!(
		#	pool_idxs_to_add, unused_full_cuts, plate_duals, threshold, n_max,
		#	sort_by_rp, rps_buffer
		#)

		pool_idxs_to_add, num_positive_rp_vars, was_above_threshold =
			_safe_recompute_idxs_to_add!(
				unused_full_cuts, plate_duals, threshold, n_max, sort_by_rp
			)
	end
	if debug
		println("qt_iters_in_iterative = $(stats.current_iter)")
		println("qt_threshold_hits_in_iterative = $(stats.qt_threshold_hits)")
	end

	return _clean_partial_byproduct!(model, lidxs, full_bp, part_bp)
end

@timeit TIMER function _final_pricing!(
	model, plate_duals, bp :: ByproductPPG2KP{D, S, P},
	LB :: Float64, LP :: Float64
) where {D, S, P}
	# TODO: check if a tolerance is needed in the comparison. Query it from the
	# solver/model if possible (or use eps?).
	kept = # the final pricing (get which variables should remain)
		@. floor(_reduced_profit(bp.cuts, (plate_duals,)) + LP) >= LB
	# yes, _delete_vars! receive which variables to keep, it is strange but it
	# will not be changed now
	new_bp = _delete_vars!(bp, model, kept, :cuts_made)
	return kept, new_bp
end

#=
function _integralize!(model, bp)
	naturally_only_binary = all(di -> di == 1, bp.d)
	extraction_vars = model[:picuts]
	if naturally_only_binary
		delete_lower_bound.(extraction_vars)
		delete_upper_bound.(extraction_vars)
		set_binary.(extraction_vars)
	else
		set_integer.(extraction_vars)
	end
	set_integer.(model[:cuts_made])

	return
end
=#


# This is a little hacky, but this option exists just to not throw away some
# code that took hard work, and to be able to include some statements in the
# paper about how growing rows does not make the code faster unfortunately.
# Ideally, this should be an option of _create_restricted_byproduct, but
# it is saner to keep it separated.
function _create_partially_restricted_byproduct(
	full_bp :: ByproductPPG2KP{D, S, P}, r_full_cut_idxs :: AbstractVector
) :: Tuple{ByproductPPG2KP{D, S, P}, BPLinkage{P}} where {D, S, P}
	# Summary: The cuts are equal to a restricted model. The plates and
	# extractions may differ: all extraction are included; if some extraction
	# is not in the restricted model, then the plate was also not in the
	# restricted model, and so both extraction and plate are included.

	# * The cuts are the `only_reachable = true` cuts. As the lidxs. No change.
	new_bp, lidxs = _create_restricted_byproduct(
		full_bp; given_rc_idxs = r_full_cut_idxs, only_reachable = true
	)
	# Get the full indexes of the extractions not in the restricted model.
	# NOTE: this is only relevant because in the revised model the extractions
	# can happen directly from plates not generated by restricted cuts.
	# In the original Furini model this is unnecessary, new_exts_idxs will
	# be empty in such case.
	new_exts_idxs = _setdiff_sorted(keys(full_bp.np), lidxs.exts_part2full)
	isempty(new_exts_idxs) && return new_bp, lidxs
	# Get the extractions (from the full_bp, with "wrong"/full plate indexes).
	new_exts = full_bp.np[new_exts_idxs]
	# Get the full plate indexes accessed from those extractions.
	plis_of_new_exts = sort!(unique!(map(first, new_exts)))
	# Get the full plate indexes not yet in in the model.
	missing_plis = setdiff(plis_of_new_exts, lidxs.plis_part2full)
	# Add the missing plates to the partial byproduct.
	append!(new_bp.pli_lwb, full_bp.pli_lwb[missing_plis])
	# Update the linked indexes, so we can use `plis_part2full` after.
	_append_to_part!(lidxs.plis_part2full, lidxs.plis_full2part, missing_plis)
	# Update the extractions to use partial plate indexes instead of full.
	map!((pli, pii) -> (lidxs.plis_full2part[pli], pii), new_exts, new_exts)
	# Add the extractions to the byproduct.
	append!(new_bp.np, new_exts)
	# Update the linked indexes of the extractions.
	_append_to_part!(lidxs.exts_part2full, lidxs.exts_full2part, new_exts_idxs)

	return new_bp, lidxs
end

# SEE SECTION 4.2 of Furini 2016: unfortunately they decided to complicate the
# method further by adding two parameters alpha and beta. The first and third
# items of the first list of section 4.4 complicate things further. In fact,
# both the greedy heuristic and the restricted model are solved but with a time
# limit (the best solution found in the middle of this process is used for the
# true final pricing).
@timeit TIMER function _furini_pricing!(
	model, full_bp, p, start :: Float64, options :: Dict{String, Any}
) where {D, S, P}
	# First let us unpack what we need from the options.
	bm :: BaseModel = (options["faithful2furini2016"] ? FURINI : BECKER)
	heuristic_seed :: Int = options["heuristic-seed"]
	verbose :: Bool = options["verbose"] & !options["quiet"]
	# both "expected" and "guaranteed" mean 'true' in this context
	mip_start :: Bool = (options["MIP-start"] != "none")
	switch_method :: Int = options["Gurobi-LP-method-inside-furini-pricing"]
	limit :: Float64 = options["building-time-limit"]
	sort_by_rp :: Bool = !options["no-sort-in-iterative-pricing"]
	alpha, beta = options["pricing-alpha"], options["pricing-beta"]
	minimal_subset = options["minimal-subset-in-iterative-pricing"]

	d = full_bp.d # shorten the name

	timings = TimeSection[]
	verbose && push!(timings, "restricted_pricing_time")
	if JuMP.solver_name(model) == "Gurobi" && switch_method > -2
		old_method = get_optimizer_attribute(model, "Method")
		# Change the LP-solving method to dual simplex, this allows for better
		# reuse of the partial solution. Could only be set after
		# `_restricted_final_pricing!` but this just cause a new set of problems.
		set_optimizer_attribute(model, "Method", switch_method)
	end

	# Builds the restricted model, with smaller number of cuts, extractions,
	# and plate constraints.
	r_full_cut_idxs = _all_restricted_cuts_idxs(full_bp)
	restricted_bp, lidxs = _create_restricted_byproduct(
		full_bp; given_rc_idxs = r_full_cut_idxs
	)
	r_inv_idxs = VarInvIndexes(restricted_bp)
	_build_base_model!(model, p, restricted_bp, r_inv_idxs, options)
	if verbose # Necessary for experiments.
		println("qt_cmvars_restricted = $(length(restricted_bp.cuts))")
		println("qt_pevars_restricted = $(length(restricted_bp.np))")
		println("qt_plates_restricted = $(length(restricted_bp.pli_lwb))")
	end

	# Use heuristic solution as LB, and relaxation as UB, to cut variables
	# from the restricted model (by a _final_pricing!), and the solve it
	# to get the LB for the unrestricted model.
	bkv, raw_ws_full, early_return, priced_r_bp = _restricted_final_pricing!(
		model, p, restricted_bp, lidxs, heuristic_seed, verbose, mip_start,
		bm, start, limit
	)
	LB = convert(Float64, bkv)

	if verbose
		println("unrestricted_LB = $LB")
		println("solved_priced_restricted_model = $(!early_return)")
		close_time = close_and_print!(timings, "restricted_pricing_time")
		push!(timings, TimeSection("iterative_pricing_time", close_time))
	end

	# If those fail, we need may need to rethink the iterative pricing,
	# because the use of max_profit seems to assume no negative profit items.
	@assert all(e -> e >= zero(e), d)
	@assert all(e -> e >= zero(e), p)

	# The iterative pricing could be better explained in the original paper.
	# Our interpretation is that it starts with all plate constraints, and all
	# extractions (`y` variables), but only with the cuts (`x` variables) that
	# represent restricted cuts. However, different from the restricted model,
	# these restricted cut variables include cuts over plates that are not
	# generated by other cuts present, in other words, cuts over 'unreachable'
	# plates are allowed. If `minimal_subset` is enabled, then only the plates
	# extracted by extractions or generated by restricted cuts are included (and,
	# consequently, only the restricted cuts over restricted-obtaineable plates
	# exist). In such case, plates/rows are, too, dinamically added to the model
	# (not only cut/columns).
	empty!(model)
	iterative_bp, lidxs = if minimal_subset
		_create_partially_restricted_byproduct(full_bp, r_full_cut_idxs)
	else
		_create_restricted_byproduct(
			full_bp; given_rc_idxs = r_full_cut_idxs, only_reachable = false
		)
	end
	full_inv_idxs = VarInvIndexes(full_bp)
	ip_inv_idxs = VarInvIndexes(iterative_bp)
	_build_base_model!(
		model, p, iterative_bp, ip_inv_idxs, options; build_LP_not_MIP = true
	)
	if verbose
		# After _iterative_pricing! these values can have changed.
		println("qt_cmvars_before_iterated = $(length(iterative_bp.cuts))")
		println("qt_pevars_before_iterated = $(length(iterative_bp.np))")
		println("qt_plates_before_iterated = $(length(iterative_bp.pli_lwb))")
	end
	final_iter_bp = _iterative_pricing!(
		model, full_bp, iterative_bp, lidxs, full_inv_idxs.pli2pair,
		sum(p .* d), alpha, beta, sort_by_rp, verbose, start, limit
	)
	#_check_linkage(lidxs)
	LP = objective_value(model)
	if verbose
		# After _iterative_pricing! these values can have changed.
		println("unrestricted_UB = $LP")
		println("qt_cmvars_after_iterated = $(length(final_iter_bp.cuts))")
		println("qt_pevars_after_iterated = $(length(final_iter_bp.np))")
		println("qt_plates_after_iterated = $(length(final_iter_bp.pli_lwb))")
		close_time = close_and_print!(timings, "iterative_pricing_time")
		push!(timings, TimeSection("final_pricing_time", close_time))
	end

	# If the upper bound (LP) computed by `_iterative_pricing!` proves that the
	# already found lower bound (LB) is optimal, then we do not need to finish
	# the model building and can just return the optimal solution. The check
	# below assumes all pieces have integer prices/'profit values' (they
	# can be Float64, they just need to not have a fractionary part for this
	# to be correct).
	@assert LB ≈ round(LB, RoundNearest)
	if (LP - LB) < 1.0
		optimum = _get_cut_pattern(raw_ws_full..., full_bp, verbose)
		return FOUND_OPTIMUM, ModelByproduct(full_bp, optimum)
	end

	# Create an array with the value of the constraint duals in the full
	# model. Throw away the old model and create the full model.
	# Use the duals to make a final pricing of the full model and finally
	# get the Priced PP-G2KP.
	full_plate_duals = zeros(Float64, length(full_bp.pli_lwb))
	full_plate_duals[lidxs.plis_part2full] .= dual.(model[:plate_cons])
	empty!(model)
	_build_base_model!(model, p, full_bp, full_inv_idxs, options)
	final_kept, final_bp = _final_pricing!(
		model, full_plate_duals, full_bp, LB, LP
	)

	if verbose
		# After _final_pricing! these values can have changed.
		println("qt_cmvars_after_final = $(length(final_bp.cuts))")
		println("qt_pevars_after_final = $(length(final_bp.np))")
		println("qt_plates_after_final = $(length(final_bp.pli_lwb))")
		close_time = close_and_print!(timings, "final_pricing_time")
		push!(timings, TimeSection("pricing_final_ws_time", close_time))
	end
	throw_if_timeout_now(start, limit)

	if mip_start
		# Shift the full indexes considering the cut variables deleted by
		# _final_pricing!, MIP-start the model with the restricted solution.
		shift_idxs!(raw_ws_full[3], final_kept)
		raw_mip_start!(model, raw_ws_full...)
	end

	# NOTE: the flag description says that Gurobi's Method is changed just for
	# _iterative_pricing! but here we change it back only after _final_pricing!,
	# the reasoning follows: _final_pricing! does not call optimize! and is not
	# affected by the LP-solving Method selected, so it is not a problem changing
	# the parameter only after it; changing the parameter before _final_pricing!
	# discards the duals of the last solve inside _iterative_pricing! and would
	# force _final_pricing! to solve the relaxation again from scratch.
	if JuMP.solver_name(model) == "Gurobi" && switch_method > -2
		# Undo the change done before (look at where old_method is defined).
		set_optimizer_attribute(model, "Method", old_method)
		#set_optimizer_attribute(model, "Presolve", -1)
	end

	verbose && close_and_print!(timings, "pricing_final_ws_time")
	return BUILT_MODEL, ModelByproduct(final_bp)
end

# TODO: for now it takes the seed and use it with the same heuristic as
# the furini pricing, but in the future we need to allow it to receive
# an arbitrary LB value, and let it call the heuristic if it is not passed
# NOTE: the p is necessary exactly because of the heuristic.
# NOTE: the _becker_pricing! can be called over both a FURINI or a BECKER
# model, this is the parameter bm.
@timeit TIMER function _becker_pricing!(
	model, bp :: ByproductPPG2KP{D, S, P}, p, start :: Float64,
	options :: Dict{String, Any}
) where {D, S, P}
	# First let us unpack what we need from the options.
	bm = (options["faithful2furini2016"] ? FURINI : BECKER) :: BaseModel
	heuristic_seed = options["heuristic-seed"]
	verbose = options["verbose"] & !options["quiet"]
	limit :: Float64 = options["building-time-limit"]
	# both "expected" and "guaranteed" mean 'true' in this context
	mip_start = (options["MIP-start"] != "none") :: Bool

	# Build the mode and relax it.
	build_complete_model(model, p, bp, start, options)
	pe_vars = model[:picuts]
	cm_vars = model[:cuts_made]
	all_vars = [pe_vars; cm_vars]
	@timeit TIMER "relax!" all_scvs = relax!(all_vars)

	# Get a lower bound (and possibly use it to MIP start the model).
	if mip_start
		(bkv, _, _), _ = mip_start_by_heuristic!(
			model, bp, p, heuristic_seed, bm
		)
	else
		bkv, _, _ = fast_iterated_greedy(
			bp.d, p, bp.l, bp.w, bp.L, bp.W, Xoroshiro128Plus(seed)
		)
	end
	LB = convert(Float64, bkv)

	# Solve the relaxation (i.e., get an upper bound).
	verbose && println("MARK_BECKER_PRICING_LP_SOLVE")
	@timeit TIMER "lp_solve" optimize_within_time_limit!(model, start, limit)
	@assert has_duals(model) # Check if the model is solved.
	LP = objective_value(model) # the model should be relaxed and solved now
	@timeit TIMER "sel_del_res" begin
		pe_rcvs = @. dual(LowerBoundRef(pe_vars))
		cm_rcvs = @. dual(LowerBoundRef(cm_vars))
		to_keep_pe_bits = @. LP - pe_rcvs >= LB
		to_keep_cm_bits = @. LP - cm_rcvs >= LB
		bp = _delete_vars!(bp, model, to_keep_pe_bits, :picuts)
		bp = _delete_vars!(bp, model, to_keep_cm_bits, :cuts_made)
		to_keep_all_bits = [to_keep_pe_bits; to_keep_cm_bits]
		@views restore!(all_vars[to_keep_all_bits], all_scvs[to_keep_all_bits])
	end
	return bp
end

