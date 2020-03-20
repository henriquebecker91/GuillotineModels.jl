module PPG2KP
# REVISED PLAN FOR IMPLEMENTING THE PRICING PART:
# 0) Implement the pricing of the revised model.
# 0.1) wrap the body of the no_arg_check_build_model inside another method
#   that just builds the Complete Model (with the reductions); call this
#   methods inside no_arg_check_build_model and then check if it will be
#   delivered this way or it will be changed by the pricing. Create a flag for
#   the initial phase of pricing and test if it reaches the 'if' correctly.
# 0.2) add documentation to the complete_model_build guaranteeing that the
#   first plate_amount constraints correspond to the plates of the model.
# 0.3) inside the pricing branch of no_arg_check_build_model, fix to zero all
#   non-restricted variables (how to do this?), relax all the other variables.
#   Solve this relaxed model. Save the UB.
# 0.4) call the heuristic and get the bkv.
# 0.6) do the final pricing step over the restricted model (i.e.,
#   the variables that could only be used to obtain a solution worse than
#   the heuristic are fixed to zero).
# 0.7) unrelax the variables that are not fixed to zero.
# 0.8) Solve the MIP of the priced restricted model, save it as the new bkv
#   (unless the solving process finished by time and the heuristic bkv is
#   better).
# DONE UNTIL THERE
# 1) Implement the initial pricing.
# 1.1) inside the pricing branch of no_arg_check_build_model, relax again
#   all the restricted cuts (the non-restriced variables should be fixed
#   to zero by now). Solve this relaxed model. Save the UB.
# 1.2) traverse all the variables and query the dual values of the plates
#   involved, use it to decide if the variable will or not be unfixed and
#   relaxed. If some variable was unfixed, repeat both previous step and this
#   step. NOTE: the check to know if the variable will be unfixed or not is in
#   p.9 first half of second column and in the first paragraph of p.13.
# 2) Implement the final pricing.
# 2.1) Use the UB and LB to do the final pricing (this maybe can reuse code
#   already developed), this will unfix some variables.
# 2.2) Finally, check which variables are fixed to zero and remove them both
#   from the model as from the ByproductPPG2KP. We need to check if we will
#   remove constraints that become irrelevant, if this is done we will need
#   to update all the cuts to refer to the new plate indexes.
# 3) Implement the warm-start using the heuristic.
# 3.1) Plan this section. Note that the warm-start of the final priced model
#   should happen inside the read_build_solve_print method and therefore it
#   needs to be a generic method to be implemented for each different model.
#   There may be previous warm starts to improve the model building process,
#   those will happen inside the build_model method.
# PLAN FOR IMPLEMENTING THE PRICING PART:
# After each step, document and test.
# 0) Remember to fix the parameters in the get_cut_pattern method,
#    it probably should be responsability of the return type to
#    store the types used, so they do not need to be in the parameters.
# 0.5) NOTE: It is never clear if the solution of the greedy heuristic is
#    used for warming up the model (in fact, this step is never described)
#    or if the LB is given to the model for cutting the search for lower
#    or equal values, or even used in any way that is not the final pricing.
# 1) Create a warm-start that is specific for the faithful2furini2016
#    and, consequently, will serve as base for the more complex one.
# 2) Create a considerably complex warm-start (comment the code
#    extensively) that uses the same heuristic to warm-start the
#    variant that is not faithful2furini2016. Check if the simpler
#    warm-start may (and should) be abandoned.
# 3) wrap the body of the no_arg_check_build_model inside another method
#    that just builds the Complete Model (with the reductions) but not
#    the priced one; call this methods inside no_arg_check_build_model and
#    then check if it will be delivered this way or it will be changed by
#    the pricing. Create a flag for the priced model and test if it reaches
#    the 'if' correctly.
# 4) Implement the iterative pricing method. The variables not in the
#    restricted model should be fixed to zero. The model must be solved
#    relaxed (check if a method is available for that, or if we will use
#    a reworked version of the handmade relaxing methods). What is said in
#    the paper kinda makes sense, but does not seem it is possible to
#    compute reduced costs for the entire pool based on the plates of the
#    restricted cut set. Write code that admits this may be impossible and
#    print the impossible situation and then abort. At the end of the process
#    (after the final pricing), really delete the not used variables, and
#    check on the ByproductPPG2KP to delete their correspondences.
#    SEE SECTION 4.2: unfortunately they decided to complicate the method
#    even more adding two parameters alpha and beta.
#    List items 1 and 3 (not numbered) of the first list of section 4.4
#    complicate things further. In fact, both the greedy heuristic and the
#    restricted model are solved but with a time limit (the best solution
#    found in the middle of this process is used for the true final pricing).

include("./Heuristic.jl")
include("./Args.jl")

include("./Enumeration.jl")
using .Enumeration
export ByproductPPG2KP # re-export ByproductPPG2KP from Enumeration

using ..Utilities
import ..CutPattern # type returned by get_cut_pattern
using RandomNumbers.Xorshifts
import ..to_pretty_str # for debug
import ..cut_pattern_profit # for debug

global_p = nothing # TODO: for debug, remove after

# TODO: Check if this method will be used, now that variable deletion is
# efficient (because our PR to JuMP).
# NOTE: by design, this method do not change the vars vector itself, but
# instead just calls delete over part of its contents (what mark such content
# as invalid)
#=
function delete_vars_by_pricing!(model, lb :: Float64)
	# All reduced costs need to be queryied and stored before we start modifying
	# the model.
	all_vars_time = @elapsed (vars = all_variables(model))
	@show all_vars_time
	obj = objective_value(model)
	rcs_time = @elapsed (rcs = reduced_cost.(vars))
	@show rcs_time
	n = length(vars)
	mask = Vector{Bool}()
	del_loop_time = @elapsed for i in 1:n
		var, rc = vars[i], rcs[i]
		@assert rc <= obj
		keep = obj - rc >= lb
		#del_time = 0.0
		!keep && #=(del_time = @elapsed=# delete(model, var)#)
		#@show del_time
		push!(mask, keep)
	end
	@show del_loop_time
	return mask
end
=#

# Execute Furini's 2016 Final Pricing. The returned value is a boolean list,
# with each index corresponding to a variable in all_variables(model),
# and with a true value if the variable should be kept, and false if it should
# be removed. The model parameter is expected to be a solved continuous
# relaxation of the problem (so the objective_value and the reduced_cost may
# be extracted); if this is not the case, pass the correct `ub` by means
# of the optional parameter.
function which_vars_to_keep(
	vars, lb :: Float64, ub :: Float64 = objective_value(model)
)
	# We keep the ones equal to lb (not strictly necessary) to guarantee the
	# warm start will work (the variables needed for it should be there).
end

using JuMP
using TimerOutputs
import MathOptInterface
const MOI = MathOptInterface

# INTERNAL METHOD USED ONLY IN get_cut_pattern
# If the "pattern" is the extraction of a single piece return either:
# (1) a CutPattern representing the piece, if it is the same size as the
# original plate; otherwise (2) a CutPattern of the original plate containing
# a CutPattern representing the piece.
function extraction_pattern(
	bmr :: ByproductPPG2KP{D, S, P}, e_idx :: Int
) :: CutPattern{D, S} where {D, S, P}
	pli, pii = bmr.np[e_idx] # the plate index and the piece index
	L, W = bmr.pli_lwb[pli] # the plate dimensions
	piece = CutPattern(bmr.l[pii], bmr.w[pii], pii)
	L == bmr.l[pii] && W == bmr.w[pii] && return piece
	return CutPattern(L, W, false, CutPattern{D, S}[piece])
end

# INTERNAL METHOD USED ONLY IN get_cut_pattern
# Given a vector of JuMP variables, return a vector of their indexes in the
# original vector and another vector storing the rounded non-zero variable
# values. (The kwargs are passed to the `isapprox` method, to define what is
# considered to be "non-zero" in a floating-point world.)
# NOTE: takes the model only to be able to check if the variable is valid.
function gather_nonzero(
	model :: JuMP.Model, vars :: Vector{JuMP.VariableRef}; kwargs...
) :: NTuple{2, Vector{Int}}
	idxs = Vector{Int}()
	vals = Vector{Int}()
	for (idx, var) in enumerate(vars)
		!is_valid(model, var) && continue
		val = value(var)
		isapprox(val, 0.0; kwargs...) && continue
		push!(idxs, idx)
		push!(vals, round(Int, val, RoundNearest))
	end
	idxs, vals
end

# INTERNAL USE.
# Check if a single piece is extracted from the original plate.
function check_if_single_piece_solution(
	np :: Vector{Tuple{P, D}}, nzpe_idxs :: Vector{Int}
) :: Int where {D, P}
	single_extraction_idx = 0
	for e_idx in nzpe_idxs
		pli = np[e_idx][1]
		!isone(pli) && continue
		single_extraction_idx = e_idx
	end

	return single_extraction_idx
end

# INTERNAL USE.
# Find the index in `cuts` of the first cut that has `pp` as its parent plate
# and has a positive `num_uses_in_sol` too; return the index after
# decrementing the `num_uses_in_sol[index]` by one. Returns `nothing` if
# there is no plate that satisfy such conditions.
function consume_cut(
  cuts :: Vector{NTuple{3, P}}, num_uses_in_sol :: Vector{D}, pp :: P
) :: Union{Nothing, P} where {D, P}
	for i in keys(cuts)
		if num_uses_in_sol[i] > zero(D) && first(cuts[i]) == pp
			num_uses_in_sol[i] -= one(D)
			return convert(P, i)
		end
	end
	return nothing
end

# INTERNAL USE.
# Given a list of the cuts used in a solution, the number of times their appear
# in the solution, and the index of the root cut in such list, return a list of
# cut indexes in topological ordering (any non-root cut over some plate only
# appears if a previous cut has made a copy of that plate type available).
function build_cut_idx_stack(
	nz_cuts :: Vector{NTuple{3, P}}, qt_cuts :: Vector{D}, root_cut_idx :: Int
) :: Vector{Int} where {D, P}
	cut_idx_stack = Vector{Int}()
	push!(cut_idx_stack, root_cut_idx)
	cuts_available = deepcopy(qt_cuts)
	next_cut = 1
	while next_cut <= length(cut_idx_stack)
		_, fc, sc = nz_cuts[cut_idx_stack[next_cut]]
		@assert !iszero(fc)
		fc_idx = consume_cut(nz_cuts, cuts_available, fc)
		fc_idx !== nothing && push!(cut_idx_stack, fc_idx)
		if !iszero(sc)
			sc_idx = consume_cut(nz_cuts, cuts_available, sc)
			sc_idx !== nothing && push!(cut_idx_stack, sc_idx)
		end
		next_cut += 1
	end

	return cut_idx_stack
end

# INTERNAL USE.
# NOTE: only patterns is modified.
# Given the non-zero (i.e., used) piece extractions (and the number of
# times they are used), add them to the patterns dictionary (that
# translates a plate id to a list of patterns for which the root node
# is a plate of that type, in this case, we are adding the patterns
# that are just piece extractions).
function add_used_extractions!(
	patterns :: Dict{Int64, Vector{CutPattern{D, S}}},
	nzpe_idxs :: Vector{Int},
	nzpe_vals :: Vector{Int},
	bmr :: ByproductPPG2KP{D, S, P}
) :: Dict{Int64, Vector{CutPattern{D, S}}} where {D, S, P}
	for (i, np_idx) in enumerate(nzpe_idxs)
		pli, pii = bmr.np[np_idx]
		#@show np_idx
		#@show pli, pii
		#@show nzpe_vals[i]
		#@show bmr.l[pii], bmr.w[pii]
		extractions = extraction_pattern.(bmr, repeat([np_idx], nzpe_vals[i]))
		if haskey(patterns, pli)
			append!(patterns[pli], extractions)
		else
			patterns[pli] = extractions
		end
	end
	patterns
end

# INTERNAL USE.
# NOTE: only patterns is modified.
# Starting from the cuts closest to the piece extractions (i.e.,
# non-leaf nodes closest to a leaf node) start building the CutPattern
# (tree) bottom-up.
function bottom_up_tree_build!(
	patterns :: Dict{Int64, Vector{CutPattern{D, S}}},
	nz_cut_idx_stack :: Vector{Int},
	nz_cuts :: Vector{NTuple{3, P}},
	nz_cuts_ori :: BitArray{1},
	bmr :: ByproductPPG2KP{D, S, P}
) :: Dict{Int64, Vector{CutPattern{D, S}}} where {D, S, P}
	for cut_idx in reverse(nz_cut_idx_stack)
		#@show cut_idx
		pp, fc, sc = nz_cuts[cut_idx]
		#@show pp, fc, sc
		child_patts = Vector{CutPattern{D, S}}()
		if !iszero(fc) && haskey(patterns, fc) && !isempty(patterns[fc])
			push!(child_patts, pop!(patterns[fc]))
			isempty(patterns[fc]) && delete!(patterns, fc)
		end
		if !iszero(sc) && haskey(patterns, sc) && !isempty(patterns[sc])
			push!(child_patts, pop!(patterns[sc]))
			isempty(patterns[sc]) && delete!(patterns, sc)
		end
		ppl, ppw = bmr.pli_lwb[pp][1], bmr.pli_lwb[pp][2]
		#@show ppl, ppw
		!haskey(patterns, pp) && (patterns[pp] = Vector{CutPattern{D, S}}())
		push!(patterns[pp], CutPattern(
			ppl, ppw, nz_cuts_ori[cut_idx], child_patts
		))
	end

	patterns
end

import ..get_cut_pattern
function get_cut_pattern(
	model_type :: Val{:PPG2KP}, model :: JuMP.Model, ::Type{D}, ::Type{S},
	build_model_return :: ByproductPPG2KP{D, S}
) :: CutPattern{D, S} where {D, S}
	# 1. Check if there can be an extraction from the original plate to a
	#    single piece. If it may and it happens, then just return this single
	#    piece solution; otherwise the first cut is in `cuts_made`, find it
	#    (i.e., traverse non-zero cuts_made variables checking for a hcut or
	#    vcut that has the original plate as the parent plate).
	# 2. Create a topological ordering of the cuts starting from the root cut.
	# 3. Initialize a pool of plate type "uses"/patterns with the piece
	#    extractions in the solution (i.e., leaf nodes).
	# 4. Traverse the reverse of the topological ordering and build the
	#    'pattern tree' in a bottom-up fashion. The children of each pattern
	#    are queryed from the pool of plate type "uses", and removed from it
	#    to be only present inside their parent pattern. At the end of the
	#    process there should remain a single pattern in the pool that is the
	#    root pattern.
	bmr = build_model_return
	pe = model[:picuts] # Piece Extractions
	cm = model[:cuts_made] # Cuts Made
	# We need to define an absolute tolerance, as solvers often return values
	# very similar to zero (when they should return zero) and this depends on
	# the instance problem coefficients.
	# MOI.get(backend(model), MOI.RelativeGap()) was initially used, but
	# it is too much of a headache, GLPK sometimes does not believe it is
	# available because we solved a linear model before.
	atol = eps(objective_value(model)) # hope this is reasonable

	# non-zero {piece extractions, cuts made} {indexes,values}
	nzpe_idxs, nzpe_vals = gather_nonzero(model, pe; atol = atol)
	nzcm_idxs, nzcm_vals = gather_nonzero(model, cm; atol = atol)
	#@show nzpe_idxs
	#@show nzpe_vals
	#@show nzcm_idxs
	#@show nzcm_vals

	sps_idx = check_if_single_piece_solution(bmr.np, nzpe_idxs)
	!iszero(sps_idx) && return extraction_pattern(bmr, sps_idx)

	# The cuts actually used in the solution.
	sel_cuts = bmr.cuts[nzcm_idxs]
	#@show sel_cuts
	# If the cut in `sel_cuts` is vertical or not.
	ori_cuts = nzcm_idxs .>= bmr.first_vertical_cut_idx
	# The index of the root cut (cut over the original plate) in `sel_cuts`.
	root_idx = findfirst(cut -> isone(cut[1]), sel_cuts)
	root_idx === nothing && return CutPattern(bmr.L, bmr.W, zero(D))

	cut_idx_stack = build_cut_idx_stack(sel_cuts, nzcm_vals, root_idx)

	# `patterns` translates a plate index (pli) into a list of all "uses" of that
	# plate type. Such "uses" are CutPattern objects, either pieces or more
	# complex patterns. When we discover a plate is used in the solution such
	# "use" is added to the known uses in `patterns`, and we arbitrarily
	# associate it to some "uses" of its children (which are then are removed
	# from `patterns`). Consequently, at the end, the `patterns` should have only
	# the key of the original plate, and its associated value should be a Vector
	# of a single CutPattern.
	patterns = Dict{Int64, Vector{CutPattern{D, S}}}()

	# Insert all piece extractions (i.e., plate to piece patterns) into
	# `patterns`.
	add_used_extractions!(patterns, nzpe_idxs, nzpe_vals, bmr)

	bottom_up_tree_build!(patterns, cut_idx_stack, sel_cuts, ori_cuts, bmr)

	global global_p
	if !isone(length(patterns)) || !haskey(patterns, 1) || !isone(length(patterns[1]))
		z = 0
		z_root = 0
		for (key, subpatts) in patterns
			@show key
			for subpatt in subpatts
				println(to_pretty_str(subpatt))
				x = cut_pattern_profit(subpatt, global_p)
				key == 1 && (z_root = x)
				z += x
			end
		end
		@show z
		@show z_root
	end

	@assert isone(length(patterns))
	@assert haskey(patterns, 1)
	@assert isone(length(patterns[1]))

	return patterns[1][1]
end

# TODO: Check why this method is used if a structure like SortedLinkedLW
# would answer this more efficiently and be aware of the type used for the
# index.
"""
    min_l_fitting_piece(l, w, L, W)

!!! **Internal use.**

Given a plate `L`x`W` and two pieces dimensions list (`l`, `w`),
return the index of the piece of smallest length that fits the
plate (the width dimension may preclude this from being just
the piece of smallest length).
"""
function min_l_fitting_piece(l, w, L, W)
	@assert length(l) == length(w)
	n = length(l)
	min_i = zero(n)
	min_l = zero(eltype(l))
	for i = one(n):n
		(l[i] > L || w[i] > W) && continue
		if iszero(min_i) || l[i] < min_l
			min_i = i
			min_l = l[i]
		end
	end
	min_i
end

"""
    search_approx_cut(pp, fcl, fcw, max_diff, approx_l, nnn, pli_lwb) :: P

!!! **Internal use.**

"""
function search_approx_cut(
	pp :: P, # parent plate
	fcl :: S, # first child length
	fcw :: S, # first child width
	max_diff :: S, # how much smaller the approx dimension may be
	approx_l :: Bool, # if true, length is approx, otherwise, w is approx
	nnn :: Vector{NTuple{3, P}},
	pli_lwb :: Vector{Tuple{S, S, P}},
) :: P where {D, S, P}
	for i in one(P):convert(P, length(nnn))
		(pp_, fc, _) = nnn[i]
		pp_ != pp && continue
		if approx_l
			pli_lwb[fc][1] <= fcl && pli_lwb[fc][1] >= (fcl - max_diff) &&
			pli_lwb[fc][2] == fcw && return i
		else
			pli_lwb[fc][2] <= fcw && pli_lwb[fc][2] >= (fcw - max_diff) &&
			pli_lwb[fc][1] == fcl && return i
		end
	end
	@show pp
	@show pli_lwb[pp]
	@show fcl
	@show fcw
	@show max_diff
	@assert false # this should not be reachable
end

"""
    search_cut_or_symmetry(pp, fcl, fcw, nnn, pli_lwb)

!!! **Internal use.**

"""
function search_cut_or_symmetry(
	pp :: P, # parent plate
	fcl :: S, # first child length
	fcw :: S, # first child width
	nnn :: Vector{NTuple{3, P}},
	pli_lwb :: Vector{Tuple{S, S, P}},
) :: P where {D, S, P}
	pll, plw, _ = pli_lwb[pp]
	@assert fcl == pll || fcw == plw
	for i in one(P):convert(P, length(nnn))
		(pp_, fc, _) = nnn[i]
		pp_ != pp && continue
		pli_lwb[fc][1] == fcl && pli_lwb[fc][2] == fcw && return i
	end
	@assert false # this should not be reachable
end

"""
    search_cut(pp, fcl, fcw, nnn, pli_lwb)

!!! **Internal use.**

"""
function search_cut(
	pp :: P, # parent plate
	fcl :: S, # first child length
	fcw :: S, # first child width
	nnn :: Vector{NTuple{3, P}},
	pli_lwb :: Vector{Tuple{S, S, P}},
) :: P where {D, S, P}
	for i in one(P):convert(P, length(nnn))
		(pp_, fc, _) = nnn[i]
		pp_ != pp && continue
		pli_lwb[fc][1] == fcl && pli_lwb[fc][2] == fcw && return i
	end
	@assert false # this should not be reachable
end

# Warm start faithful2furini.
function warm_start_f2f(
	model, l, w, L, W,
	pat :: Vector{Vector{D}},
	bp :: ByproductPPG2KP{D, S, P}
	#round2disc wait to see if this is needed
	# which other model building options will need to be passed to this?
) where {D, S, P}

end

# TODO: check if the piece will be immediatelly extracted or does it need
# an extra cut to make a plate small enough to extract the piece?
# TODO: if the piece needs an extra trim, the things are even worse than
# I thought at first, there is no guarantee the trim will be obtained
# in a single cut. It may be necessary to make many fake trim cuts to
# get a plate of the right size to make an extraction. It is just much
# simpler to add fake extraction variables for every plate that would need
# fake trims than to add all cuts needed to arrive at it. On the other side,
# this is something interesting to comment on the paper, the lack of guarantee
# on the number of cuts needed to extract the pieces from the plates, a
# 'downside' that seems not to be a problem for the solver, of have common
# worst-cases.
# TODO: comment all code this method is specially hard to follow
# TODO: the warm-start for the flag faithful2furini enabled and disabled will
# need to be different? the rules for which plates exist and which do not are
# different from one to another.
# NOTE: this method only work for simple patterns in which:
# (1) the cuts are two-staged (i.e., the pattern is justa a vector of vector);
# (1.5) by consequence, the cuts are restricted;
# (2) the first stage cuts vertically (width strips);
# (3) the first piece of each strip gives the width of the strip;
# If you need to warm-start with a more complex pattern, create another
# method with the same name, and another type for parameter `pat`.
"""
    warm_start(model, l, w, L, W, pat, pli_lwb, nnn, np; faithful2furini = false)

!!! **Internal use.**
"""
function warm_start(
	model, l, w, L, W,
	pat :: Vector{Vector{D}},
	pli_lwb :: Vector{Tuple{S, S, P}},
	nnn :: Vector{NTuple{3, P}},
	np :: Vector{Tuple{P, D}};
	faithful2furini2016 = false
	#round2disc wait to see if this is needed
	# which other model building options will need to be passed to this?
) where {D, S, P}
	@assert !faithful2furini2016
	# the initial residual plate is L, W
	# visit the outer vector in reverse
	# if the current head of stripe is smaller than half residual plate
	# then search for a vertical cut on the residual plate, with the right width
	#   for the first child and enable it, change the residual plate to
	#   be the second child
	# else assert this is the last stripe, just use the remaining plate (second
	#   child of the last cut, or the whole root if there is just one stripe)
	# after finishing the strip processing, for each plate that is a stripe:
	#   do the same as the first stage, but for the subplate and the opposite cut
	#   orientation (as the inner vectors are not sorted by length, we need to
	#   sort them ourselves);
	# finally, for every subplate that will become a piece, connect it to a piece
	#   for faithful2furini2016 we need to trim the plate and have it with exact
	#   plate size; for !faithful2furini2016 we just link directly to np
	rw = W # remaining width, initialized with root plate width
	rpli = 1 # NOTE: the root plate is guaranteed to have index 1
	cut_var_vals = Dict{P, Int}() # the amount of times each cut was made
	first_stage_plates = Vector{P}() # the plate index of all stripes
	rpat = reverse(pat)
	final = false
	for stripe in rpat
		@assert !final
		@assert !isempty(stripe)
		ws = w[first(stripe)]
		if ws <= div(rw, 2)
			cut = search_cut(rpli, L, ws, nnn, pli_lwb)
			cut_var_vals[cut] = 1 + get(cut_var_vals, cut, 0)
			_, fc, rpli = nnn[cut]

			rw -= ws
			push!(first_stage_plates, fc)
			fcl, fcw = pli_lwb[fc]
			println("0\ti\t$(fcl)\t$(fcw)")
		else
			# If the stripe width is already more than half residual plate, then it
			# is the last stripe (the stripes are ordered in increase-or-keep order,
			# and cannot exist a larger width stripe in the remaining space).
			push!(first_stage_plates, rpli)
			scl, scw = pli_lwb[rpli]
			println("0\tl\t$(scl)\t$(scw)")
			final = true
		end
	end
	# TODO: both loops need to be by index, and stop before the last element,
	# the last element need to be checked against the min_piece_l/w and the
	# remaining space in that dimension, there may be needed to make a cut
	# that creates an unused first child and a used second child.
	# The elements before the last do not have this problem because the
	# ordering guarantee that they are smaller than half plate and therefore
	# a cut for them exist.
	# TODO: Consider: would it be simpler to just create missing variables?
	# this would work for any model type, and would make this code unaffedcted
	# by new variable reductions; the increase in the number of variables
	# would be very small.

	println("start of second stage")
	pli2pii_qt = Dict{Tuple{P, D}, Int}()
	for i in 1:length(first_stage_plates)
		# For each stripe, reset the remaining info.
		rpli = first_stage_plates[i] # get the index of the plate/stripe
		rl = pli_lwb[rpli][1] # get the length of the plate/stripe
		ws = pli_lwb[rpli][2] # stripe width (do not diminish)

		pli2pii = Vector{Tuple{P, D}}()
		for piece in sort(rpat[i], by = pii -> l[pii])
			lp = l[piece]
			wp = w[piece]
			if lp <= div(rl, 2)
				cut = search_cut(rpli, lp, ws, nnn, pli_lwb)
				cut_var_vals[cut] = 1 + get!(cut_var_vals, cut, 0)
				_, fc, rpli = nnn[cut]

				fcl, fcw, _ = pli_lwb[fc]
				println("$(i)\ti\t$(fcl)\t$(fcw)")
				println("$(piece)\tp\t$(lp)\t$(wp)")

				min_w = min_l_fitting_piece(w, l, fcw, fcl)
				if faithful2furini2016
				elseif fcw >= wp + min_w
					# The piece may have a small width relative to the plate width
					# (i.e., the piece of smallest width could fit there) if this is
					# the case, we cannot extract it directly from tne plate.
					# There are two options: (1) if the piece is smaller than half
					# the plate, it can be cut as first child; (2) otherwise,
					# we need to make an extra cut (a fake trim) in which
					# the second child will be a plate "not too big" to allow a
					# direct extraction of the piece from it.
					if wp <= div(fcw, 2)
						# Cut the fc again, at the piece width this time, use the
						# first child to extract the piece.
						cut2 = search_cut(fc, lp, wp, nnn, pli_lwb)
						cut_var_vals[cut2] = 1 + get!(cut_var_vals, cut2, 0)
						_, fc2, _ = nnn[cut2]

						fc2l, fc2w, _ = pli_lwb[fc2]
						println("$(i)\ti\t$(scl)\t$(scw)")
						println("$(piece)\tp\t$(lp)\t$(wp)")

						pli2pii_qt[(fc2, piece)] = 1 + get(pli2pii_qt, (fc2, piece), 0)
					else
						trimw = fcw - wp
						# There is no guarantee of existence of a cut in which the
						# second child will have the exact size of the piece, however
						# there is guarantee that a cut between the exact size and
						# the exact size - (min_w - 1) exists. Note that cuts smaller
						# than this lower bound do not interest us, as they would have
						# the same problem again (the second child would need a fake
						# trim cut again).
						cut2 = search_approx_cut(
							rpli, lp, trimw, min_w - 1, false, nnn, pli_lwb
						)
						cut_var_vals[cut2] = 1 + get!(cut_var_vals, cut2, 0)
						_, _, sc = nnn[cut2]

						scl, scw, _ = pli_lwb[sc]
						println("$(i)\ti\t$(scl)\t$(scw)")
						println("$(piece)\tp\t$(lp)\t$(wp)")

						pli2pii_qt[(sc, piece)] = 1 + get(pli2pii_qt, (sc, piece), 0)
					end
				else
					pli2pii_qt[(fc, piece)] = 1 + get(pli2pii_qt, (fc, piece), 0)
				end

				rl -= lp
			else # if it is the last piece of a strip
				pll, plw, _ = pli_lwb[rpli]
				@assert plw == ws
				min_l_pii = min_l_fitting_piece(l, w, pll, plw)
				min_l = iszero(min_l_pii) ? pll : l[min_l_pii]
				if pll >= lp + min_l
					triml = pll - lp
					cut = search_approx_cut(
						rpli, triml, ws, min_l - 1, true, nnn, pli_lwb
					)
					cut_var_vals[cut] = 1 + get!(cut_var_vals, cut, 0)
					_, _, sc = nnn[cut]

					scl, scw, _ = pli_lwb[sc]
					println("$(i)\tl\t$(scl)\t$(scw)")
					println("$(piece)\tp\t$(lp)\t$(wp)")
					rpli = sc
				end
				pll, plw, _ = pli_lwb[rpli]
				@assert plw == ws
				min_w_pii = min_l_fitting_piece(w, l, plw, pll)
				min_w = iszero(min_w_pii) ? plw : l[min_w_pii]
				if plw >= wp + min_w
					trimw = plw - wp
					cut = search_approx_cut(
						rpli, pll, trimw, min_w - 1, false, nnn, pli_lwb
					)
					cut_var_vals[cut] = 1 + get!(cut_var_vals, cut, 0)
					_, _, sc = nnn[cut]

					scl, scw, _ = pli_lwb[sc]
					println("$(i)\tl\t$(scl)\t$(scw)")
					println("$(piece)\tp\t$(lp)\t$(wp)")
					rpli = sc
				end

				pli2pii_qt[(rpli, piece)] = 1 + get(pli2pii_qt, (rpli, piece), 0)
			end
		end
	end
	for (var_index, var_value) in cut_var_vals
		set_start_value(model[:cuts_made][var_index], var_value)
	end
	if faithful2furini2016
		# HERE WE NEED TO DO THE TRIM CUT IF FAITHFUL2FURINI IS ENABLED
	else
		for (pair, value) in pli2pii_qt
			@show pli_lwb[pair[1]]
			@show l[pair[2]], w[pair[2]]
			picut_idx = findfirst(isequal(pair), np)
			@show picut_idx
			set_start_value(model[:picuts][picut_idx], value)
		end
	end

	model
end

# HIGH LEVEL EXPLANATION OF THE MODEL
#
# Variables:
#
# `picuts[n, pii]`: Integer. The number of pieces `pii` generated from
#   subplates of type `n`.
# `cuts_made[n1, n2, n3]`: Integer. The number of subplates of type
#   `n1` that are cut into subplates `n2` and `n3` (horizontal and
#   vertical cuts are together for now).
#
# Objective function:
#
# Maximize the profit of the pieces cut.
#   sum(p[pii] * picuts[_, pii])
#
# Constraints:
#
# There is exactly one of the original plate, which may be used for cutting
# or extracting a piece.
#   sum(picuts[1, _]) + sum(cuts_made[1, _,  _]) <= 1
# The number of subplates available depends on the number of plates that have
# it as children.
#   sum(picuts[n1>1, _]) + sum(cuts_made[n1>1, _, _]) <=
#     sum(cuts_made[_, n2, n3])
#     where n2 == n1 or n3 == n1, doubling cuts_made[_, n2, n3] if n2 == n3
# The number of pieces of some type is always less than or equal to the demand.
#   sum(picuts[_, pii]) <= d[pii]
#
# Unnecessary constraints:
#
# The number of times a pair pli-pii appear is at most the min between:
# d[pii] and the number of subplates pli that fit in the original plate.
#   sum(picuts[n, pii]) <= min(d[pii], max_fits[n])
"""
    build_complete_model(model, d, p, l, w, L, W; options)

TODO: describe method.
Note: this method guarantee that the first M constraints of type
`(GenericAffExpr{Float64,VariableRef}, LessThan{Float64})` in the model,
where M is the number of plates, will correspond to the plates of the
model (i.e., as by the field ByproductPPG2KP.pli_lwb).
"""
function build_complete_model(
	model, d :: Vector{D}, p :: Vector{P}, l :: Vector{S}, w :: Vector{S},
	L :: S, W :: S, options :: Dict{String, Any} = Dict{String, Any}()
) :: ByproductPPG2KP{D, S, P} where {D, S, P}
	@timeit "enumeration" begin
	@assert length(d) == length(l) && length(l) == length(w)
	num_piece_types = convert(D, length(d))

	sllw = SortedLinkedLW(D, l, w)
	# TODO: change enumeration to also use a Dict?
	byproduct = gen_cuts(P, d, sllw, L, W;
		ignore_2th_dim = options["ignore-2th-dim"],
		ignore_d = options["ignore-d"],
		round2disc = options["round2disc"],
		no_cut_position = options["no-cut-position"],
		no_redundant_cut = options["no-redundant-cut"],
		no_furini_symmbreak = options["no-furini-symmbreak"],
		faithful2furini2016 = options["faithful2furini2016"]
	)
	hvcuts, pli_lwb, np = byproduct.cuts, byproduct.pli_lwb, byproduct.np
	num_plate_types = length(pli_lwb)

	# Order the plate-piece pairs by my ranking of importance: how much
	# absolute are is wasted. Other rankings include: how much relative area
	# is wasted; how much is the profit of the a squared unit of the plate used;
	# group them by pii just for making it easier to the demand constraint.
	#=sort!(np, lt = function((pli1, pii1), (pli2, pii2))
		pli1l, pli1w, _ = pli_lwb[pli1]
		pli2l, pli2w, _ = pli_lwb[pli2]
		pii1l, pii1w = l[pii1], w[pii1]
		pii2l, pii2w = l[pii2], w[pii2]
		convert(P, pli1l - pii1l) * (pli1w - pii1w) < (
		convert(P, pli2l - pii2l) * (pli2w - pii2w))
	end)=#

	# pli2pair: inverse index which given a plate index will return a list of
	# all picuts indexes (np variable) that cut some piece from some plate.
	pli2pair = [Vector{P}() for _ = 1:num_plate_types]
	# pii2pair: the same as pli2pair but for piece indexes.
	pii2pair = [Vector{P}() for _ = 1:num_piece_types]

	# The vectors below all have the same length as the number of plate types (to
	# allow indexing by plate type). The value of a position is a vector of
	# arbitrary length and irrelevant index, the values of this inner vector are
	# cut indexes. Such cut indexes are related to the plate that is the index of
	# the outer vector.
	# any CHILD plate to respective CUT indexes
	child2cut = [Vector{P}() for _ = 1:num_plate_types]
	# any PARENT plate to respective CUT indexes
	parent2cut = [Vector{P}() for _ = 1:num_plate_types]

	# Initialize all inverse indexes.
	for i in eachindex(hvcuts)
		parent, fchild, schild = hvcuts[i]
		push!(parent2cut[parent], i)
		push!(child2cut[fchild], i)
		#=if iszero(schild)
			@show parent
			@show pli_lwb[parent]
			@show fchild
			@show pli_lwb[fchild]
			@show schild
		end=#
		#@assert faithful2furini2016 || !iszero(schild)
		!iszero(schild) && push!(child2cut[schild], i)
	end

	for i in eachindex(np)
		pli, pii = np[i]
		push!(pli2pair[pli], i)
		push!(pii2pair[pii], i)
	end
	end # @timeit "enumeration_related", old: time_to_enumerate_plates

	@timeit "JuMP_calls" begin
	# If all pieces have demand one, a binary variable will suffice to make the
	# connection between a piece type and the plate it is extracted from.
	naturally_only_binary = all(di -> di <= 1, d)
	if naturally_only_binary
		@variable(model, picuts[1:length(np)], Bin)
	else
		#@variable(model, picuts[1:length(np)] >= 0, Int)
		@variable(model,
			0 <= picuts[i = 1:length(np)] <= min(pli_lwb[np[i][1]][3], d[np[i][2]]),
		Int)
	end

	@variable(model, cuts_made[1:length(hvcuts)] >= 0, Int)

	# The objective function is to maximize the profit made by extracting
	# pieces from subplates.
	@objective(model, Max,
		sum(p[pii] * sum(picuts[pii2pair[pii]]) for pii = 1:num_piece_types)
	)

	# c1: There is just one of the original plate, and so it can be only used
	# to extract a single piece xor make a single cut that would make two new
	# subplates available.
	@constraint(model,
		sum(picuts[pli2pair[1]]) + sum(cuts_made[parent2cut[1]]) <= 1
	)

	# c2: for each subplate type that is not the original plate, such subplate
	# type will be available the number of times it was the child of a cut,
	# subtracted the number of times it had a piece extracted or used for
	# further cutting.
	for pli in 2:num_plate_types
		@constraint(model,
			sum(picuts[pli2pair[pli]]) + sum(cuts_made[parent2cut[pli]]) <=
			sum(cuts_made[child2cut[pli]])
		)
	end

	if options["use-c25"]
		# c2.5: The amount of each subplate type generated by cuts (and used either
		# as a piece or as a intermediary plate) is bounded by the amount that can be
		# cut from the original plate.
		for pli in 2:num_plate_types
			@constraint(model,
				sum(picuts[pli2pair[pli]]) + sum(cuts_made[parent2cut[pli]]) <=
				pli_lwb[pli][3]
			)
		end
	end

	# c3: the amount of each piece type extracted from different plate types
	# cannot surpass the demand for that piece type.
	for pii in 1:num_piece_types
		@constraint(model, sum(picuts[pii2pair[pii]]) <= d[pii])
	end

	lb = options["lower-bound"]
	if !iszero(lb)
		@constraint(model,
			sum(p[pii]*sum(picuts[pii2pair[pii]]) for pii = 1:num_piece_types) >= (lb + 1)
		)
	end

	ub = options["upper-bound"]
	if ub < sum(d .* p)
		@constraint(model,
			sum(p[pii]*sum(picuts[pii2pair[pii]]) for pii = 1:num_piece_types) <= ub
		)
	end
	end # @timeit "calls_to_JuMP", old: time_to_solver_build

	return byproduct
end

function all_restricted_cuts_idxs(
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
function setdiff_sorted(u, s)
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
function reduced_profit(
	cut :: Tuple{P, P, P}, constraints :: Vector{T}
) :: Float64 where {P, T}
	pp, fc, sc = cut
	sc_dual = iszero(sc) ? 0.0 : dual(constraints[sc])
	# The reduced profit computation shown in the paper is done here.
	return dual(constraints[pp]) - dual(constraints[fc]) - sc_dual
end

function restricted_final_pricing(
	model, rc_idxs :: Vector{Int}, d :: Vector{D}, p :: Vector{P},
	bp :: ByproductPPG2KP{D, S, P}, rng, debug :: Bool = false
) #=:: Tuple{P, CutPattern{D, S}, Vector{SavedVarConf}} =# where {D, S, P}
	# First we save all the variables bounds and types and relax all of them.
	n = num_variables(model)
	pe = model[:picuts] # Piece Extractions
	cm = model[:cuts_made] # Cuts Made
	debug && begin
		@show n
		@show length(pe)
		@show length(cm)
	end
	# This assert is here because, if this becomes false some day, the code
	# will need to be reworked. Now it only relax/fix variables in those sets
	# so it will need to be sensibly extended to new sets of variables.
	@assert n == length(cm) + length(pe)
	pe_svcs = save_and_relax!.(pe)
	cm_svcs = save_and_relax!.(cm)
	#@assert length(cm) == length(bp.cuts)
	rc_vars = cm[rc_idxs]
	# Even if they are contained in all_svcs, save the var config of the
	# restricted cuts again, because we will need to restore only them after
	# and finding them in all_svcs is a drag.
	rc_svcs = SavedVarConf.(rc_vars)
	# Every variable that is not a restricted cut is a unrestricted cut.
	uc_idxs = setdiff_sorted(keys(cm), rc_idxs)
	# Fix all unrestricted cuts (pool vars).
	uc_vars = cm[uc_idxs]
	fix.(uc_vars, 0.0; force = true)
	# Solve the relaxed restricted model.
	@timeit "relaxed_restricted" optimize!(model)
	#@show JuMP.dual_status(model)
	#@show JuMP.primal_status(model)
	# Get the solution of the heuristic.
	bkv, _, shelves = Heuristic.iterated_greedy(
		d, p, bp.l, bp.w, bp.L, bp.W, rng
	)
	sol = Heuristic.shelves2cutpattern(shelves, bp.l, bp.w, bp.L, bp.W)
	LB = convert(Float64, bkv)
	#@show LB
	# Check if everything seems ok with the values obtained.
	if termination_status(model) == MOI.OPTIMAL
		#@show objective_bound(model)
		#@show objective_value(model)
		UB = objective_value(model)
		#@show value.(rc_vars)
		!isinf(UB) && (UB = floor(UB + eps(UB))) # theoretically sane UB
		# Note: the iterated_greedy solve the shelf version of the problem,
		# that is restricted, so it cannot return a better known value than the
		# restricted upper bound.
		@assert LB <= UB
	else
		@error(
			"For some reason, solving the relaxed restricted model did not" *
			" terminate with optimal status (the status was" *
			" $(termination_status(model))). The code is not prepared to deal" *
			" with this possibility and will abort."
		)
		exit(1)
	end
	rc_cuts = bp.cuts[rc_idxs]
	plate_cons = all_constraints(
		model, GenericAffExpr{Float64, VariableRef}, MOI.LessThan{Float64}
	)
	# TODO: check if a tolerance is needed in the comparison. Query it from the
	# solver/model if possible (or use eps?).
	# Do the final pricing, get which variables should be removed.
	rc_fix_bitstr = @. floor(UB - reduced_profit(rc_cuts, (plate_cons,))) <= LB
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
	fix.(rc_vars[rc_fix_bitstr], 0.0; force = true)
	# The discrete variables are the restricted cuts that were not removed
	# by the final pricing of the restricted model and will be in the MIP
	# of the restricted model.
	rc_discrete_bitstr = .!rc_fix_bitstr
	discrete_vars = rc_vars[rc_discrete_bitstr]
	restore!.(discrete_vars, rc_svcs[rc_discrete_bitstr])
	# All piece extractions also need to be restored.
	restore!.(pe, pe_svcs)
	@timeit "discrete_restricted" optimize!(model) # restricted MIP solved
	# If some primal solution was obtained, compare its value with the
	# heuristic and keep the best one (the model can give a worse value
	# as any variables that cannot contribute to a better solution are
	# disabled, and the heuristic may be already optimal).
	if primal_status(model) == MOI.FEASIBLE_POINT
		model_obj = objective_value(model)
		model_obj = floor(model_obj + eps(model_obj))
		#@show model_obj
		if model_obj > LB
			bkv = convert(P, model_obj)
			sol = get_cut_pattern(Val(:PPG2KP), model, D, S, bp)
		end
	end
	# The variables are left all relaxed but with only the variables used
	# in the restricted priced model unfixed, the rest are fixed to zero.
	relax!.(rc_vars[rc_discrete_bitstr])
	relax!.(pe)

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
# Return: the number of positive reduced profit variables found.
function recompute_idxs_to_add!(
	pool_idxs_to_add, pool, plate_cons, threshold, n_max :: P,
	debug :: Bool = false
) :: P where {P}
	found_above_threshold = false
	num_positive_rp_vars = zero(P)
	# TODO: consider if the optimization of using the pool_idxs_to_add as a
	# buffer and deciding between overwrite and push! is worth.
	empty!(pool_idxs_to_add)

	for (pool_idx, cut) in pairs(pool)
		rp = reduced_profit(cut, plate_cons)
		rp > 0.0 && continue # non-positive reduced profit is irrelevant to us
		num_positive_rp_vars += 1
		if rp > threshold
			!found_above_threshold && empty!(pool_idxs_to_add)
			found_above_threshold = true
			push!(pool_idxs_to_add, pool_idx)
			# If n_max variables above the threshold exist, only them are used.
			length(pool_idxs_to_add) >= n_max && break
		elseif !found_above_threshold
			# Unfortunately we cannot stop here if we find n_max variables because
			# we can find one above the threshold later yet.
			push!(pool_idxs_to_add, pool_idx)
		end
	end

	debug && @show found_above_threshold
	# Only n_max variables are added/unfixed in a single iteration.
	n_max <= length(pool_idxs_to_add) && resize!(pool_idxs_to_add, n_max)

	return num_positive_rp_vars
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
function iterative_pricing(
	model, cuts :: Vector{NTuple{3, P}}, #bp :: ByproductPPG2KP{D, S, P},
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

	# We use copy, and not deepcopy, because we do not want the original
	# containers to change, the contents either need to change or cannot change.
	unused_cuts = copy(cuts)
	unfixed_at_start = (!is_fixed).(model[:cuts_made])
	unused_vars = copy(model[:cuts_made])
	# We just need the unused vars/cuts in the pricing process, once a variable
	# is "added" (in truth, unfixed) it is always kept, so it does not need to be
	# reevaluated.
	deleteat!(unused_vars, unfixed_at_start)
	deleteat!(unused_cuts, unfixed_at_start)
	@assert length(unused_cuts) == length(unused_vars)
	# The p̄ value in the paper, if there are variables with reduced profit
	# above this threshold then they are added (and none below the threshold),
	# but the n_max limit of variables added in a single iteration is respected.
  threshold = max_profit * beta
	@assert threshold > zero(threshold) # threshold is always larger than zero
	# The vectors are all allocated here (but used inside `price!`) for
	# reuse (avoiding reallocating them every loop).
	to_unfix = Vector{eltype(keys(cuts))}()
	# The first 'number of plates' constraints here are the constraints that
	# relate to respective 1:'number of plates' plates.
	plate_cons = all_constraints(
		model, GenericAffExpr{Float64, VariableRef}, MOI.LessThan{Float64}
	)
	optimize!(model) # the last solve before this was MIP and has no duals
	# Do the initial pricing, necessary to compute n_max, and that is always done
	# (i.e., the end condition can only be tested after this first loop).
	num_positive_rp_vars = recompute_idxs_to_add!(
		to_unfix, unused_cuts, plate_cons, threshold, typemax(P), debug
	)
	# The maximum number of variables unfixed at each iteration.
	n_max = round(P, num_positive_rp_vars * alpha, RoundUp)
	debug && @show num_positive_rp_vars
	debug && @show n_max
	# TODO: CONTINUE HERE EXPLAIN RESIZE BELOW
	resize!(to_unfix, n_max)
	# The iterative pricing continue until there are no variables
	# with positive reduced profit value to unfix.
	while !iszero(num_positive_rp_vars)
		unfix.(unused_vars[to_unfix])
		deleteat!(unused_vars, to_unfix)
		deleteat!(unused_cuts, to_unfix)
		@assert length(unused_cuts) == length(unused_vars)
		optimize!(model)
		num_positive_rp_vars = recompute_idxs_to_add!(
			to_unfix, unused_cuts, plate_cons, threshold, n_max, debug
		)
		debug && @show num_positive_rp_vars
		debug && @show length(to_unfix)
	end
	return
end

function final_pricing(
	model, byproduct :: ByproductPPG2KP{D, S, P},
	options :: Dict{String, Any} = Dict{String, Any}()
) :: ByproductPPG2KP{D, S, P} where {D, S, P}
	return byproduct
end

function Base.empty!(model::Model)::Model
    # The method changes the Model object to, basically, the state it was when
    # created (if the optimizer was already pre-configured). The two exceptions
    # are: optimize_hook and operator_counter. One is left alone because it is
    # related to optimizer attributes and the other is just a counter for a
    # single-time warning message (so keeping it helps to discover
    # inneficiencies).
    MOI.empty!(model.moi_backend)
    empty!(model.shapes)
    model.nlp_data = nothing
    empty!(model.obj_dict)
    empty!(model.ext)
    return model
end

function no_arg_check_build_model(
	model, d :: Vector{D}, p :: Vector{P}, l :: Vector{S}, w :: Vector{S},
	L :: S, W :: S, options :: Dict{String, Any} = Dict{String, Any}()
) :: ByproductPPG2KP{D, S, P} where {D, S, P}
	global global_p = p
	byproduct = build_complete_model(model, d, p, l, w, L, W, options)
	options["no-pricing"] && return byproduct

	rng = Xoroshiro128Plus(options["pricing-heuristic-seed"])
	rc_idxs = all_restricted_cuts_idxs(byproduct)
	UB, sol, pe_svcs, cm_svcs = restricted_final_pricing(
		model, rc_idxs, d, p, byproduct, rng, options["debug"]
	)
	# If those fail, we need may need to rethink the iterative pricing,
	# because the use of max_profit seems to assume no negative profit items.
	@assert all(e -> e >= zero(e), d)
	@assert all(e -> e >= zero(e), p)
	# TODO: check if is needed to pass the whole bp structure
	@timeit "iterative_pricing" iterative_pricing(
		model, byproduct.cuts, sum(d .* p),
		options["pricing-alpha"], options["pricing-beta"], options["debug"]
	)
	LB = objective_value(model)

	# TODO: the empty! and the uncommented return are just while the pricing
	# is not implemented.
	empty!(model)
	return build_complete_model(model, d, p, l, w, L, W, options)
	#return byproduct
end

#=
# This block of code needs to find a home here, as it is model specific.
# Before, it inhabited the script code.
	if p_args["warm-start"]
		@assert !p_args["faithful2furini2016"] && !p_args["flow-model"]
		(heur_obj, heur_sel, heur_pat), heur_time, _, _, _ = @timed iterated_greedy(
			d, p, l, w, L, W, MersenneTwister(p_args["solver-seed"])
		)
		if !no_csv_out
			@show heur_time
			@show heur_obj
			@show heur_sel
			@show heur_pat
		end
		saved_bounds = PPG2KP.disable_unrestricted_cuts!(
			m, sort(l), sort(w), hvcuts, pli_lwb
		)
		heur_ws_time = @elapsed PPG2KP.warm_start(
			m, l, w, L, W, heur_pat, pli_lwb, hvcuts, np;
			faithful2furini2016 = p_args["faithful2furini2016"]
		)
		!no_csv_out && @show heur_ws_time
		restricted_time = @elapsed optimize!(m)
		!no_csv_out && @show restricted_time
		restricted_sol = value.(all_variables(m))
		restricted_objval = objective_value(m)
		PPG2KP.restore_bound!.(saved_bounds)
		set_start_value.(all_variables(m), restricted_sol)
	end

	# HERE IT WAS THE DIVISION OF BEFORE BUILDING THE MODEL AND AFTER BUILDING IT

	vars_before_deletes = all_variables(m)
	if p_args["final-pricing"] || p_args["relax2lp"]
		original_settings = relax_all_vars!(m)
	end
	time_to_solve_model = @elapsed optimize!(m)
	if (p_args["final-pricing"] || p_args["relax2lp"]) && !no_csv_out
		println("time_to_solve_relaxed_model = $(time_to_solve_model)")
	end
	if p_args["final-pricing"]
		# get the given lower bound on the instance if there is one
		best_lb = get(p_args["lower-bounds"], instfname_idx, 0)
		# if warm-start is enabled and the lb is better, use it
		p_args["warm-start"] && (best_lb = max(best_lb, restricted_objval))
		# delete all variables which would only reduce obj to below the lower bound
		#if !no_csv_out
		#  @profile (mask = delete_vars_by_pricing!(m, Float64(best_lb)))
		#  Profile.print()
		#else
		#  mask = delete_vars_by_pricing!(m, Float64(best_lb))
		#end
		kept = which_vars_to_keep(m, Float64(best_lb))
		unrelax_vars!(vars_before_deletes, m, original_settings)
		saved_bounds = fix_vars!(all_variables(m)[.!kept])
		if !no_csv_out
			num_vars_after_pricing = sum(kept)
			num_vars_del_by_pricing = length(kept) - num_vars_after_pricing
			@show num_vars_after_pricing
			@show num_vars_del_by_pricing
		end
		time_to_solve_model += @elapsed optimize!(m)
	end
	sleep(0.01) # shamefully needed to avoid out of order messages from cplex
=#

import ..Utilities.Args.create_normalized_arg_subset
import ..Utilities.Args.accepted_arg_list
import ..build_model
"""
    build_model(::Val{:PPG2KP}, model, d, p, l, w, L, W[, options])


"""
function build_model(
	::Val{:PPG2KP}, model, d :: Vector{D}, p :: Vector{P},
	l :: Vector{S}, w :: Vector{S}, L :: S, W :: S,
	options :: Dict{String, Any} = Dict{String, Any}()
) :: ByproductPPG2KP{D, S, P} where {D, S, P}
	norm_options = create_normalized_arg_subset(
		options, accepted_arg_list(Val(:PPG2KP))
	)
	return no_arg_check_build_model(model, d, p, l, w, L, W, norm_options)
end

end # module

