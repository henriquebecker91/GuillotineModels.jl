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
# an iterator over all the elements in `u` but not in `s`.
function _setdiff_sorted(u, s)
	isempty(s) && return deepcopy(u)
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

#=
struct BPLinkage{P}
	part2full
	full2part
end
=#

@timeit TIMER function _create_restricted_byproduct(
   bp :: ByproductPPG2KP{D, S, P}
) :: ByproductPPG2KP{D, S, P} where {D, S, P}
	# The cuts/plates/extractions of a restricted model are obtained by a
	# two-step process: (1) get only the cuts in which the cut length/width
	# match the corresponding dimension of a piece; (2) detect which plates,
	# extractions, and restricted cuts are reachable from those cuts.
	# Note that the cuts obtained in (1) can occur in plates that cannot be
	# obtained from the cuts in (1), consequently, the number of restricted cuts
	# often drastically reduces in (2) (and the values become similar to the ones
	# obtained by Furini).
	fvc = if bp.first_vertical_cut_idx > length(bp.cuts)
		nothing
	else
		bp.cuts[bp.first_vertical_cut_idx]
	end
	rc_idxs = _all_restricted_cuts_idxs(bp)
	rc_fvci = searchsortedfirst(rc_idxs, bp.first_vertical_cut_idx)
	rc = bp.cuts[rc_idxs]
	qt_re, re_bits, qt_rc, rc_bits, qt_rp, rp_bits = _reachable(
		bp.np, rc, lastindex(bp.pli_lwb)
	)
	re = deleteat!(copy(bp.np), .!re_bits)
	rp = deleteat!(copy(bp.pli_lwb), .!rp_bits)
	deleteat!(rc, .!rc_bits) # rc is already a copy of a subset
	rc_fvci = rc_fvci > length(rc_idxs) ? length(rc) + 1 : sum(rc_bits[1:rc_fvci])
	# If this assert fails, it means our index tricks using searchsortedfirst
	# and sum have failed us.
	@assert rc_fvci > length(rc) || rc[rc_fvci] == fvc
	# After the assert, we can update the values of the cuts and extractions to
	# refer to the new plates.
	qt_old_pl = length(bp.pli_lwb)
	oldpli2newpli = zeros(qt_old_pl)
	setindex!.((oldpli2newpli,), 1:length(rp), (1:qt_old_pl)[rp_bits])
	# TODO: put these loops inside a internal method of their own, I have
	# the feeling I have already written them many times.
	for (i, (pp, fc, sc)) in pairs(rc)
		new_pp, new_fc = oldpli2newpli[pp], oldpli2newpli[fc]
		new_sc = iszero(sc) ? sc : oldpli2newpli[sc]
		@assert new_pp > 0 && new_pp <= length(rp)
		if !(new_fc > 0 && new_fc <= length(rp))
			@show i
			@show pp
			@show fc
			@show sc
			@show oldpli2newpli[fc]
			@show length(rp)
		end
		@assert new_sc >= 0 && new_fc <= length(rp)
		rc[i] = (new_pp, new_fc, new_sc)
	end
	for (i, (plate, piece)) in pairs(re)
		@assert !iszero(oldpli2newpli[plate])
		re[i] = (oldpli2newpli[plate], piece)
		@assert re[i][1] > 0 && re[i][1] <= length(rp)
	end

	# FUTURE PLAN: if we start returning the bitmasks to keep the link
	# between original and restricted byproducts, then we can use
	# broadcasted getindex and setindex! to create a bitmask for cuts using
	# rc_idxs and rc_bits.

	return ByproductPPG2KP{D, S, P}(rc, rc_fvci, re, rp, bp.l, bp.w, bp.L, bp.W)
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
	model, d :: Vector{D}, p :: Vector{P},
	bp :: ByproductPPG2KP{D, S, P}, seed, debug :: Bool, mip_start :: Bool,
	bm :: BaseModel, start :: Float64 = time(),
	limit :: Float64 = float(60*60*24*365), options :: Dict{String, Any} = Dict{String, Any}()
) where {D, S, P}
	empty!(model)
	bp = _create_restricted_byproduct(bp)
	r_inv_idxs = VarInvIndexes(bp)
	_build_base_model!(model, d, p, bp, r_inv_idxs, options)
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

	# Debug necessary for experiments.
	if debug
		println("qt_cmvars_restricted = $(length(bp.cuts))")
		println("qt_pevars_restricted = $(length(bp.np))")
		println("qt_plates_restricted = $(length(bp.pli_lwb))")
	end

	# The 'best known value'/'primal bound' of the heuristic is needed to do
	# the pricing. If we can MIP-start the model, we can just call
	# `mip_start_by_heuristic!` and get the `bkv` returned, otherwise we need
	# to call just the heuristic to get the bkv (and throw away the rest).
	if mip_start
		(bkv, _, _), heuristic_raw_ws = mip_start_by_heuristic!(
			model, bp, d, p, seed, bm
		)
	else
		bkv :: P, _, _ = fast_iterated_greedy(
			d, p, bp.l, bp.w, bp.L, bp.W, Xoroshiro128Plus(seed)
		)
	end
	restricted_LB = LB = convert(Float64, bkv)
	debug && @show restricted_LB

	# Check if the heuristic did not blow the time limit, but only after
	# the information output that happens after it.
	throw_if_timeout_now(start, limit)

	# Solve the relaxed restricted model.
	debug && println("MARK_FURINI_PRICING_RESTRICTED_LP_SOLVE")
	@timeit TIMER "lp_solve" optimize_within_time_limit!(model, start, limit)
	# Check if everything seems ok with the values obtained.
	if termination_status(model) == MOI.OPTIMAL
		restricted_LP = LP = objective_value(model)
		# Note: the iterated_greedy solve the shelf version of the problem,
		# that is restricted, so it cannot return a better known value than the
		# restricted "upper bound"/"linear(ized) problem".
		debug && @show restricted_LP
		@assert restricted_LB <= restricted_LP
	else
		error(
			"For some reason, solving the relaxed restricted model did not" *
			" terminate with optimal status (the status was" *
			" $(termination_status(model))). The code is not prepared to deal" *
			" with this possibility and will abort."
		)
	end

	# If the assert below fails, then the solver used reaches an optimal point of
	# a LP but do not has duals for it, what should not be possible, some problem
	# has happened (for an example, the solver does not recognize the model as a
	# LP but think it is a MIP for example, CPLEX sometimes does this).
	@assert has_duals(model)

	# Now let us do the restricted pricing, using the heuristic as LB and
	# the continuous relaxation as UB.
	plate_cons = model[:plate_cons]
	# TODO: check if a tolerance is needed in the comparison. Query it from the
	# solver/model if possible (or use eps?).
	rc_discrete_bitstr = # the final pricing (get which variables should remain)
		@. floor(_reduced_profit(bp.cuts, (dual(plate_cons),)) + LP) >= LB
	# ... and which ones should be fixed to zero.
	rc_fix_bitstr = .!rc_discrete_bitstr
	if debug
		qt_cmvars_priced_restricted = sum(rc_discrete_bitstr)
		qt_cmvars_deleted_by_heuristic_pricing = length(rc_fix_bitstr) -
			qt_cmvars_priced_restricted
		@show qt_cmvars_priced_restricted
		@show qt_cmvars_deleted_by_heuristic_pricing
	end
	# TODO: here we need to be able to keep the linkage between the old BP
	# and the new one
	bp = _delete_vars!(bp, model, rc_discrete_bitstr, :cuts_made)
	restore!.(model[:picuts], pe_svcs)
	restore!.(model[:cuts_made], cm_svcs[rc_discrete_bitstr])
	# restricted MIP solved
	debug && println("MARK_FURINI_PRICING_RESTRICTED_MIP_SOLVE")
	@timeit TIMER "mip_solve" optimize_within_time_limit!(model, start, limit)
	# If we MIP start the restricted model with a feasible solution, it should
	# be impossible to get a different status here.
	@assert !mip_start || primal_status(model) == MOI.FEASIBLE_POINT
	# The assert above is valid, but this does not mean we actually deal
	# with the possibility of not having a primal solution below.
	model_obj = round(P, objective_value(model), RoundNearest)
	if debug
		restricted_stop_reason = termination_status(model)
		restricted_stop_code = Int(restricted_stop_reason)
		println("restricted_obj_value = $model_obj")
		println("restricted_obj_bound = $(objective_bound(model))")
		@show restricted_stop_reason
		@show restricted_stop_code
	end
	if mip_start && model_obj > bkv
		# Clean the old MIP start (that comes from the heuristic) and replace it
		# with the solution of the restricted model (yes, take the solution and
		# re-input it as warm-start just to have an extra guarantee it will be
		# used). The ideal would be saving it and passing it along to just be
		# applied before returning the model. This because some solvers can discard
		# the start value of all variables if some variables are deleted (even if
		# the deleted ones are not in the mip start). However, this would need us
		# to either: save JuMP variables instead of indexes; or update the indexes
		# every time a variable is deleted. For now this is too much effort
		# because Gurobi log shows that the MIP start is being recognized even with
		# variables being deleted after.
		# `shift_idxs!` below correct the indexes considering the deleted vars.
		shift_idxs!(heuristic_raw_ws[3], rc_discrete_bitstr)
		unset_mip_start!(model, heuristic_raw_ws[1], heuristic_raw_ws[3])
		raw_mip_start!(model, save_mip_start(model)...)
	end
	model_obj > bkv && (bkv = model_obj)

	# The variables are left all relaxed but with only the variables used
	# in the restricted priced model unfixed, the rest are fixed to zero.
	@timeit TIMER "relax_ss_rc" relax!(cm)
	@timeit TIMER "relax_pe" relax!(pe)

	error(
		"Commits done inside branch furinifix have no guarantee of" *
		" being functional. This version aborts at the end of the" *
		" restricted furini pricing (if it is called)."
	)

	return bkv, pe_svcs, cm_svcs
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
	pool_idxs_to_add, pool, plate_cons, threshold, n_max :: P,
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
		rps_buffer[pool_idx] = _reduced_profit(cut, plate_cons)
		# non-positive reduced profit is irrelevant to us
		rps_buffer[pool_idx] <= 0.0 && continue
		num_positive_rp_vars += one(P)
		if rps_buffer[pool_idx] > threshold
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
	model, cuts :: Vector{NTuple{3, P}}, cm_svcs :: Vector{SavedVarConf},
	max_profit :: P, alpha :: Float64, beta :: Float64, sort_by_rp :: Bool,
	debug :: Bool = false, start :: Float64 = time(),
	limit :: Float64 = float(60*60*24*365)
) :: Nothing where {P}
	# Summary of the method:
	# The variables themselves are not iterated, bp.cuts is iterated. The indexes
	# of the plates in each cut element are the indexes of the associated
	# constraints. Query the duals of the three associated constraints, compute
	# the reduced profit. Then, we need the list of idxs of variables with
	# positive reduced profit but below p̄ threshold, and the one above p̄
	# threshold. If this second list is non-empty, then we take the first n_max
	# items and unfix the respective variables. If it is empty, we do this but to
	# the first list instead. If the first list is empty too, then we finished
	# the iterative pricing process.

	# This code is simplified by the assumption all variables had lower and
	# upper bounds in their original discrete (integer or binary) configuration.
	# If the asserts below triggers, then this method need to be reevaluated.
	@assert all(svc -> svc.had_ub, cm_svcs)
	@assert all(svc -> svc.had_lb, cm_svcs)
	unused_lbs = getfield.(cm_svcs, :lb)
	unused_ubs = getfield.(cm_svcs, :ub)
	# We use copy, and not deepcopy, because we want to avoid only the original
	# *container* to change, the contents either need to change or cannot be
	# changed (i.e., are immutable).
	unused_cuts = copy(cuts)
	unfixed_at_start = (!is_fixed).(model[:cuts_made])
	if debug
		println("qt_cmvars_before_iterated = $(sum(unfixed_at_start))")
	end
	unused_vars = copy(model[:cuts_made])
	# We just need the unused vars/cuts in the pricing process, once a variable
	# is "added" (in truth, unfixed) it is always kept, so it does not need to be
	# reevaluated.
	deleteat!(unused_vars, unfixed_at_start)
	deleteat!(unused_cuts, unfixed_at_start)
	deleteat!(unused_lbs, unfixed_at_start)
	deleteat!(unused_ubs, unfixed_at_start)
	allsame(x) = all(y -> y == x[1], x)
	@assert allsame(length.((unused_vars, unused_cuts, unused_lbs, unused_ubs)))
	# The p̄ value in the paper, if there are variables with reduced profit
	# above this threshold then they are added (and none below the threshold),
	# but the n_max limit of variables added in a single iteration is respected.
  threshold = max_profit * beta
	@assert threshold > zero(threshold) # threshold is always larger than zero
	# The vectors are all allocated here (but used inside `price!`) for
	# reuse (avoiding reallocating them every loop).
	to_unfix = Vector{eltype(keys(cuts))}()
	plate_cons = model[:plate_cons]
	size_var_pool_before_iterative = length(unused_vars)
	debug && @show size_var_pool_before_iterative
	flush_all_output()
	# the last solve before this was MIP and has no duals
	debug && println("MARK_FURINI_PRICING_ITERATED_LP_SOLVE_0")
	#@timeit TIMER "solve_lp" optimize_within_time_limit!(model, start, limit)
	@timeit TIMER "solve_lp" optimize!(model)
	#v_attr_pstart = Gurobi.VariableAttribute("PStart")
	#c_attr_dstart = Gurobi.ConstraintAttribute("DStart")
	#all_vars = all_variables(model)
	#relevant_constraints = vcat(plate_cons, model[:demand_con])
	#pstart = MOI.get.((backend(model),), (v_attr_pstart,), index.(all_vars))
	#dstart = MOI.get.((backend(model),), (c_attr_dstart,), index.(relevant_constraints))
	#pstart = value.(all_vars)
	#dstart = value.(relevant_constraints)
	#vector_summary(pstart)
	#vector_summary(dstart)
	#plate_cons_dual = dual.(plate_cons)
	#println("plate_con_duals stats")
	#vector_summary(plate_cons_dual)
	flush_all_output()
	# Do the initial pricing, necessary to compute n_max, and that is always done
	# (i.e., the end condition can only be tested after this first loop).
	rps_buffer = Vector{Float64}(undef, length(unused_cuts))
	initial_num_positive_rp_vars, was_above_threshold = _recompute_idxs_to_add!(
		to_unfix, unused_cuts, plate_cons, threshold, typemax(P), sort_by_rp,
		rps_buffer
	)
	pricing_threshold_hits = was_above_threshold ? 1 : 0
	debug && @show initial_num_positive_rp_vars
	# The maximum number of variables unfixed at each iteration.
	n_max = round(P, initial_num_positive_rp_vars * alpha, RoundUp)
	debug && @show n_max
	# As we called `num_positive_rp_vars!` with `typemax(P)` instead `n_max`
	# to be able to compute `n_max` in the first place, we need to resize
	# this first list to `n_max` (in the loop below `_recompute_idxs_to_add!`
	# will do it for us).
	n_max < length(to_unfix) && resize!(to_unfix, n_max)
	# Not sure if restarting the model use the old values as start values.
	#all_vars = all_variables(model)
	num_positive_rp_vars = initial_num_positive_rp_vars
	# The iterative pricing continue until there are variables to unfix.
	qt_vars_added_back_to_LP_by_iterated = zero(P)
	qt_iters = zero(P)
	while !isempty(to_unfix)
		if debug
			qt_iters += one(P)
			println("qt_positive_rp_vars_iter_$(qt_iters) = $(num_positive_rp_vars)")
			println("above_threshold_iter_$(qt_iters) = $(was_above_threshold)")
			println("qt_vars_added_to_LP_iter_$(qt_iters) = $(length(to_unfix))")
			qt_vars_added_back_to_LP_by_iterated += length(to_unfix)
		end
		@timeit TIMER "update_bounds" begin
			unfixed_vars = unused_vars[to_unfix]
			unfix.(unfixed_vars)
			set_lower_bound.(unfixed_vars, unused_lbs[to_unfix])
			set_upper_bound.(unfixed_vars, unused_ubs[to_unfix])
		end
		@timeit TIMER "deleteat!" begin
			deleteat!(unused_vars, to_unfix)
			deleteat!(unused_cuts, to_unfix)
			deleteat!(unused_lbs, to_unfix)
			deleteat!(unused_ubs, to_unfix)
		end
		@assert allsame(length.((unused_vars, unused_cuts, unused_lbs, unused_ubs)))
		#set_start_value.(all_vars, value.(all_vars))
		flush_all_output()
		debug && println("MARK_FURINI_PRICING_ITERATED_LP_SOLVE_$(qt_iters)")
		#println("USING Gurobi.reset_model!(backend(model).inner)")
		#Gurobi.reset_model!(backend(model).inner)
		#@timeit TIMER "solve_lp" optimize_within_time_limit!(model, start, limit)
		#MOI.set.((backend(model),), (v_attr_pstart,), index.(all_vars), pstart)
		#MOI.set.((backend(model),), (c_attr_dstart,), index.(relevant_constraints), dstart)
		@timeit TIMER "solve_lp" optimize!(model)
		#pstart .= value.(all_vars)
		#dstart .= value.(relevant_constraints)
		#pstart .= MOI.get.((backend(model),), (v_attr_pstart,), all_vars)
		#dstart .= MOI.get.((backend(model),), (c_attr_dstart,), relevant_constraints)
		flush_all_output()
		num_positive_rp_vars, was_above_threshold = _recompute_idxs_to_add!(
			to_unfix, unused_cuts, plate_cons, threshold, n_max, sort_by_rp,
			rps_buffer
		)
		was_above_threshold && (pricing_threshold_hits += 1)
	end
	if debug
		println("qt_iterations_of_iterated_pricing = $qt_iters")
		@show pricing_threshold_hits
		@show qt_vars_added_back_to_LP_by_iterated
	end

	return
end

@timeit TIMER function _final_pricing!(
	model, bp :: ByproductPPG2KP{D, S, P}, LB :: Float64, LP :: Float64,
	debug :: Bool = false
) where {D, S, P}
	plate_cons = model[:plate_cons]
	vars = model[:cuts_made]
	debug && @show LP
	debug && @show LB
	# If this assert fails is because some change was done to a solved model
	# and its state was reset to unoptimized.
	@assert has_duals(model)
	plate_duals = dual.(plate_cons)
	rps = _reduced_profit.(bp.cuts, (plate_duals,))
	#=if debug
		println("stats on reduced profit values of final pricing")
		vector_summary(rps) # defined in GuillotineModels.Utilities
	end=#
	to_keep_bits = @. floor(rps + LP) >= LB # the final pricing
	return _delete_vars!(bp, model, to_keep_bits, :cuts_made), to_keep_bits
end

# SEE SECTION 4.2 of Furini 2016: unfortunately they decided to complicate the
# method further by adding two parameters alpha and beta. The first and third
# items of the first list of section 4.4 complicate things further. In fact,
# both the greedy heuristic and the restricted model are solved but with a time
# limit (the best solution found in the middle of this process is used for the
# true final pricing).
@timeit TIMER function _furini_pricing!(
	model, byproduct, d, p, seed, alpha, beta, debug, mip_start :: Bool,
	bm :: BaseModel, switch_method :: Int, sort_by_rp :: Bool,
	start :: Float64 = time(), limit :: Float64 = float(60*60*24*365),
	options :: Dict{String, Any} = Dict{String, Any}()
) where {D, S, P}
	if JuMP.solver_name(model) == "Gurobi" && switch_method > -2
		old_method = get_optimizer_attribute(model, "Method")
		# Change the LP-solving method to dual simplex, this allows for better
		# reuse of the partial solution. Could only be set after
		# `_restricted_final_pricing!` but this just cause a new set of problems.
		set_optimizer_attribute(model, "Method", switch_method)
	end
	# The two possible mip-starts occur inside `_restricted_final_pricing!`.
	bkv, pe_svcs, cm_svcs = _restricted_final_pricing!(
		model, d, p, byproduct, seed, debug, mip_start, bm, start, limit, options
	)
	debug && print_past_section_seconds(TIMER, "_restricted_final_pricing!")
	LB = convert(Float64, bkv)
	# If those fail, we need may need to rethink the iterative pricing,
	# because the use of max_profit seems to assume no negative profit items.
	@assert all(e -> e >= zero(e), d)
	@assert all(e -> e >= zero(e), p)

	_iterative_pricing!(
		model, byproduct.cuts, cm_svcs, sum(d .* p), alpha, beta, sort_by_rp,
		debug, start, limit
	)
	debug && print_past_section_seconds(TIMER, "_iterative_pricing!")
	LP = objective_value(model)
	byproduct, to_keep_bits = _final_pricing!(
		model, byproduct, LB, LP, debug
	)
	debug && print_past_section_seconds(TIMER, "_final_pricing!")
	throw_if_timeout_now(start, limit)
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
	# Restore the piece extractions and the kept variables to their original
	# configuration. Note that cuts_made is now a subset of what it was before.
	@timeit TIMER "last_restore" @views restore!(
		[model[:picuts]; model[:cuts_made]], [pe_svcs; cm_svcs[to_keep_bits]]
	)

	return byproduct
end

# TODO: for now it takes the seed and use it with the same heuristic as
# the furini pricing, but in the future we need to allow it to receive
# an arbitrary LB value, and let it call the heuristic if it is not passed
# NOTE: the d and p are necessary exactly because of the heuristic.
# NOTE: the _becker_pricing! can be called over both a FURINI or a BECKER
# model, this is the parameter bm.
@timeit TIMER function _becker_pricing!(
	model, bp :: ByproductPPG2KP{D, S, P}, d, p, seed,
	debug :: Bool, mip_start :: Bool, bm :: BaseModel,
	start :: Float64 = time(), limit :: Float64 = float(60*60*24*365)
) where {D, S, P}
	pe_vars = model[:picuts]
	cm_vars = model[:cuts_made]
	all_vars = [pe_vars; cm_vars]
	@timeit TIMER "relax!" all_scvs = relax!(all_vars)
	if mip_start
		(bkv, _, _), _ = mip_start_by_heuristic!(
			model, bp, d, p, seed, bm
		)
	else
		bkv, _, _ = fast_iterated_greedy(
			d, p, bp.l, bp.w, bp.L, bp.W, Xoroshiro128Plus(seed)
		)
	end
	LB = convert(Float64, bkv)
	debug && println("MARK_BECKER_PRICING_LP_SOLVE")
	@timeit TIMER "lp_solve" optimize_within_time_limit!(model, start, limit)
	# Check if the model is solved.
	@assert has_duals(model)
	# Check if the variables are relaxed.
	#@assert all(v -> !is_integer(v) & !is_binary(v), all_vars)
	# Assert no variable considered is fixed but have lower bounds.
	#@assert all(v -> has_lower_bound(v) & !is_fixed(v), cm_vars)
	LP = objective_value(model) # the model should be relaxed and solved now
	LB = 0.0
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
