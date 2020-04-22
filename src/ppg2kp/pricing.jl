# Gives the indexes of all elements of bp.cuts that refer to a restricted cut.
# NOTE: no extraction in bp.np can be seen as an unrestricted cut.
function _all_restricted_cuts_idxs(
	bp :: ByproductPPG2KP{D, S, P}
) :: Vector{Int} where {D, S, P}
	rc_idxs = Vector{Int}()
	usl = unique!(sort(bp.l))
	usw = unique!(sort(bp.w))
	for (idx, (_, fc, _)) in pairs(bp.cuts)
		if idx < bp.first_vertical_cut_idx
			!isempty(searchsorted(usl, bp.pli_lwb[fc][1])) && push!(rc_idxs, idx)
		else
			!isempty(searchsorted(usw, bp.pli_lwb[fc][2])) && push!(rc_idxs, idx)
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
	sc_dual = iszero(sc) ? 0.0 : dual(constraints[sc])
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

# TODO: the call to the heuristic and solving the restricted MIP model
# probably should not be here, but on no_arg_check_build_model.
@timeit TIMER function _restricted_final_pricing!(
	model, rc_idxs :: Vector{Int}, d :: Vector{D}, p :: Vector{P},
	bp :: ByproductPPG2KP{D, S, P}, rng, debug :: Bool = false
) #=:: Tuple{P, CutPattern{D, S}, Vector{SavedVarConf}} =# where {D, S, P}
	# First we save all the variables bounds and types and relax all of them.
	#all_vars = all_variables(model)
	#n = length(all_vars)
	n = num_variables(model)
	pe = model[:picuts] # Piece Extractions
	cm = model[:cuts_made] # Cuts Made
	# This assert is here because, if this becomes false some day, the code
	# will need to be reworked. Now it only relax/fix variables in those sets
	# so it will need to be sensibly extended to new sets of variables.
	@assert n == length(cm) + length(pe)
	#@timeit TIMER "relax_cm"
	cm_svcs = relax!(cm)
	#@timeit TIMER "relax_pe"
	pe_svcs = relax!(pe)
	#@assert length(cm) == length(bp.cuts)
	rc_vars = cm[rc_idxs]
	rc_svcs = cm_svcs[rc_idxs]
	# Every variable that is not a restricted cut is a unrestricted cut.
	uc_idxs = _setdiff_sorted(keys(cm), rc_idxs)
	# Fix all unrestricted cuts (pool vars).
	uc_vars = cm[uc_idxs]
	#@timeit TIMER "fix_uc"
	fix.(uc_vars, 0.0; force = true)
	# Get the solution of the heuristic.
	bkv, _, shelves = @timeit TIMER "primal_heuristic" Heuristic.iterated_greedy(
		d, p, bp.l, bp.w, bp.L, bp.W, rng
	)
	sol = Heuristic.shelves2cutpattern(shelves, bp.l, bp.w, bp.L, bp.W)
	LB = convert(Float64, bkv)
	# Solve the relaxed restricted model.
	flush_all_output()
	@timeit TIMER "lp_solve" optimize!(model)
	flush_all_output()
	# Check if everything seems ok with the values obtained.
	if termination_status(model) == MOI.OPTIMAL
		#@show objective_bound(model)
		#@show objective_value(model)
		LP = UB = objective_value(model)
		#@show value.(rc_vars)
		!isinf(UB) && (UB = floor(UB + eps(UB))) # theoretically sane UB
		# Note: the iterated_greedy solve the shelf version of the problem,
		# that is restricted, so it cannot return a better known value than the
		# restricted upper bound.
		restricted_LB = LB
		restricted_LP = LP
		debug && begin
			@show restricted_LB
			@show restricted_LP
		end
		@assert restricted_LB <= restricted_LP
	else
		@error(
			"For some reason, solving the relaxed restricted model did not" *
			" terminate with optimal status (the status was" *
			" $(termination_status(model))). The code is not prepared to deal" *
			" with this possibility and will abort."
		)
		exit(1)
	end
	# TODO: if the assert below fails, then the solver used reaches an optimal
	# point of a LP but do not has duals for it, what should not be possible,
	# some problem has happened (for an example, the solver does not recognize
	# the model as a LP but think it is a MIP for example, CPLEX does this).
	@assert has_duals(model)
	#all_vars_values = value.(all_vars) # needs to be done before any changes
	rc_cuts = bp.cuts[rc_idxs]
	plate_cons = model[:plate_cons]
	if debug
		println("stats on reduced profit values of restricted final pricing")
		rps = _reduced_profit.(rc_cuts, (plate_cons,))
		vector_summary(rps) # defined in GuillotineModels.Utilities
	end
	# TODO: check if a tolerance is needed in the comparison. Query it from the
	# solver/model if possible (or use eps?).
	rc_discrete_bitstr = # the final pricing (get which variables should remain)
		@. floor(_reduced_profit(rc_cuts, (plate_cons,)) + LP) >= LB
	rc_fix_bitstr = .!rc_discrete_bitstr
	debug && begin
		restricted_vars_removed = sum(rc_fix_bitstr)
		restricted_vars_remaining = length(rc_cuts) - restricted_vars_removed
		@show restricted_vars_removed
		@show restricted_vars_remaining
		unrestricted_vars_fixed = length(uc_idxs)
		@show unrestricted_vars_fixed
	end
	# Below the SavedVarConf is not stored because they will be kept fixed
	# for now, and when restored, they will be restored to their original
	# state (rc_svcs) not this intermediary one.
	#@timeit TIMER "fix_rc_subset"
	fix.(rc_vars[rc_fix_bitstr], 0.0; force = true)
	# The discrete variables are the restricted cuts that were not removed
	# by the final pricing of the restricted model and will be in the MIP
	# of the restricted model.
	discrete_vars = rc_vars[rc_discrete_bitstr]
	#@timeit TIMER "restore_ss_cm"
	restore!(discrete_vars, rc_svcs[rc_discrete_bitstr])
	# All piece extractions also need to be restored.
	#@timeit TIMER "restore_pe"
	restore!(pe, pe_svcs)
	#set_start_value.(all_vars, all_vars_values)
	flush_all_output()
	@timeit TIMER "mip_solve" optimize!(model) # restricted MIP solved
	flush_all_output()
	# If some primal solution was obtained, compare its value with the
	# heuristic and keep the best one (the model can give a worse value
	# as any variables that cannot contribute to a better solution are
	# disabled, and the heuristic may be already optimal).
	if primal_status(model) == MOI.FEASIBLE_POINT
		model_obj = objective_value(model)
		model_obj = round(model_obj, RoundNearest)
		#@show model_obj
		if model_obj > LB
			bkv = convert(P, model_obj)
			sol = get_cut_pattern(Val(:PPG2KP), model, D, S, bp)
		end
	end
	# The variables are left all relaxed but with only the variables used
	# in the restricted priced model unfixed, the rest are fixed to zero.
	#@timeit TIMER "relax_ss_rc"
	relax!(rc_vars[rc_discrete_bitstr])
	#@timeit TIMER "relax_pe"
	relax!(pe)
	#@assert all(v -> !is_integer(v) && !is_binary(v), all_variables(model))
	#set_start_value.(all_vars, all_vars_values) # using the values of the LP

	return bkv, sol, pe_svcs, cm_svcs
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
) :: Tuple{P, Bool} where {P}
	found_above_threshold = false
	num_positive_rp_vars = zero(P)
	# TODO: consider if the optimization of using the pool_idxs_to_add as a
	# buffer and deciding between overwrite and push! is worth.
	empty!(pool_idxs_to_add)

	for (pool_idx, cut) in pairs(pool)
		rp = _reduced_profit(cut, plate_cons)
		rp <= 0.0 && continue # non-positive reduced profit is irrelevant to us
		num_positive_rp_vars += one(P)
		if rp > threshold
			!found_above_threshold && empty!(pool_idxs_to_add)
			found_above_threshold = true
			push!(pool_idxs_to_add, pool_idx)
			# If n_max variables above the threshold exist, only them are used.
			length(pool_idxs_to_add) >= n_max && break
		elseif !found_above_threshold
			# Unfortunately we cannot stop here if we find n_max variables because
			# we can find one above the threshold later yet. At least, we can stop
			# pushing new variables to the pool.
			length(pool_idxs_to_add) < n_max && push!(pool_idxs_to_add, pool_idx)
		end
	end

	# Only n_max variables are added/unfixed in a single iteration.
	#n_max <= length(pool_idxs_to_add) && resize!(pool_idxs_to_add, n_max)

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
	max_profit :: P, alpha :: Float64, beta :: Float64, debug :: Bool = false
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
	@timeit TIMER "solve_lp" optimize!(model)
	flush_all_output()
	# Do the initial pricing, necessary to compute n_max, and that is always done
	# (i.e., the end condition can only be tested after this first loop).
	initial_num_positive_rp_vars, was_above_threshold = _recompute_idxs_to_add!(
		to_unfix, unused_cuts, plate_cons, threshold, typemax(P)
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
		debug && begin
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
		@timeit TIMER "solve_lp" optimize!(model)
		flush_all_output()
		num_positive_rp_vars, was_above_threshold = _recompute_idxs_to_add!(
			to_unfix, unused_cuts, plate_cons, threshold, n_max
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

# TODO: check if the kept variables are chosen from the unfixed variables
# of from all variables.
@timeit TIMER function _final_pricing!(
	model, bp :: ByproductPPG2KP{D, S, P}, LB :: Float64, LP :: Float64,
	debug :: Bool = false
) where {D, S, P}
	plate_cons = model[:plate_cons]
	vars = model[:cuts_made]
	#unfixed_bits = !is_fixed(vars)
	#unfixed_vars, fixed_vars = _partition_by_bits(unfixed_bits, vars)
	#unfixed_cuts = bp.cuts[unfixed_bits]
	debug && @show LP
	debug && @show LB
	# Many optional small adjusts may clean the model duals, if they are not
	# available we re-solve the LP one more time to make them available again.
	!has_duals(model) && optimize!(model)
	if debug
		println("stats on reduced profit values of final pricing")
		rps = _reduced_profit.(bp.cuts, (plate_cons,))
		vector_summary(rps) # defined in GuillotineModels.Utilities
	end
	to_keep_bits = # the final pricing (get which variables should remain)
		@. floor(_reduced_profit(bp.cuts, (plate_cons,)) + LP) >= LB
	return _delete_vars!(bp, model, to_keep_bits, :cuts_made), to_keep_bits
end

# SEE SECTION 4.2 of Furini 2016: unfortunately they decided to complicate the
# method further by adding two parameters alpha and beta. The first and third
# items of the first list of section 4.4 complicate things further. In fact,
# both the greedy heuristic and the restricted model are solved but with a time
# limit (the best solution found in the middle of this process is used for the
# true final pricing).
@timeit TIMER function _furini_pricing!(
	model, byproduct, d, p, seed, alpha, beta, debug
) where {D, S, P}
	if JuMP.solver_name(model) == "Gurobi"
		old_method = get_optimizer_attribute(model, "Method")
		# Change the LP-solving method to dual simplex, this allows for better
		# reuse of the partial solution.
		set_optimizer_attribute(model, "Method", 1)
		#set_optimizer_attribute(model, "Presolve", 0)
	end

	rng = Xoroshiro128Plus(seed)
	rc_idxs = _all_restricted_cuts_idxs(byproduct)
	bkv, sol, pe_svcs, cm_svcs = _restricted_final_pricing!(
		model, rc_idxs, d, p, byproduct, rng, debug
	)
	LB = convert(Float64, bkv)
	# If those fail, we need may need to rethink the iterative pricing,
	# because the use of max_profit seems to assume no negative profit items.
	@assert all(e -> e >= zero(e), d)
	@assert all(e -> e >= zero(e), p)
	# TODO: check if is needed to pass the whole bp structure
	_iterative_pricing!(
		model, byproduct.cuts, cm_svcs, sum(d .* p), alpha, beta, debug
	)
	LP = objective_value(model)
	# TODO: change the name of the methods to include the exclamation mark
	byproduct, to_keep_bits = _final_pricing!(
		model, byproduct, LB, LP, debug
	)
	# Restore the piece extractions and the kept variables to their original
	# configuration. Note that cuts_made is now a subset of what it was before.
	@timeit TIMER "last_restore" @views restore!(
		[model[:picuts]; model[:cuts_made]], [pe_svcs; cm_svcs[to_keep_bits]]
	)
	if JuMP.solver_name(model) == "Gurobi"
		# Undo the change done before (look at where old_method is defined).
		set_optimizer_attribute(model, "Method", old_method)
		#set_optimizer_attribute(model, "Presolve", -1)
	end

	return byproduct
end

# TODO: for now it takes the seed and use it with the same heuristic as
# the furini pricing, but in the future we need to allow it to receive
# an arbitrary LB value, and let it call the heuristic if it is not passed
# NOTE: the d and p are necessary exactly because of the heuristic
@timeit TIMER function _becker_pricing!(
	bp, model, d, p, seed, debug
)
	pe_vars = model[:picuts]
	cm_vars = model[:cuts_made]
	all_vars = [pe_vars; cm_vars]
	@timeit TIMER "relax!" all_scvs = relax!(all_vars)
	flush_all_output()
	@timeit TIMER "lp_solve" optimize!(model)
	flush_all_output()
	# Check if the variables are relaxed.
	@assert all(v -> !is_integer(v) & !is_binary(v), all_vars)
	# Check if the model is solved.
	@assert has_duals(model)
	# Assert no variable considered is fixed but have lower bounds.
	@assert all(v -> has_lower_bound(v) & !is_fixed(v), cm_vars)
	LP = objective_value(model) # the model should be relaxed and solved now
	LB = 0.0
	@timeit TIMER "primal_heuristic" begin
		rng = Xoroshiro128Plus(seed)
		bkv, _, _ = Heuristic.iterated_greedy(
			d, p, bp.l, bp.w, bp.L, bp.W, rng
		)
		LB = convert(Float64, bkv)
	end
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
