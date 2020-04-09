module PPG2KP
# 2.2) Finally, check which variables are fixed to zero and remove them both
#   from the model as from the ByproductPPG2KP. We need to check if we will
#   remove constraints that become irrelevant, if this is done we will need
#   to update all the cuts to refer to the new plate indexes.
# PLAN FOR IMPLEMENTING THE WARM-START
# 0) Adapt the warm-start flag to be PPG2KP-specific and pricing-specific
#    it will just define if we will make use of the primal solutions already
#    available to warm-start the model.
# 1) The only warm-start that really needs dedicated code is the one made from
#    the heuristic to a restricted model. The other is just setting the
#    variable MIP start to the value obtained solving the restricted model.
# 1) Create a warm-start procedure that is specific for the faithful2furini2016
#    and, consequently, will serve as base for the more complex one.
# 2) Create a considerably complex warm-start (comment the code
#    extensively) that uses the same heuristic to warm-start the
#    variant that is not faithful2furini2016. Check if the simpler
#    warm-start may (and should) be abandoned.

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

using JuMP
using TimerOutputs
import MathOptInterface
const MOI = MathOptInterface

# total_number_of_plate_types should be greater-than-or-equal-to the largest
# value found in any field of cuts.
function _reachable_plate_types(
	cuts :: Vector{NTuple{3, P}}, total_number_of_plate_types :: P
) where {P}
	children = [P[] for _ in 1:total_number_of_plate_types]
	for (pp, fc, sc) in cuts
		push!(children[pp], fc)
		!iszero(sc) && push!(children[pp], sc)
	end
	reached = falses(total_number_of_plate_types)
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

#= disabled because it is O(n^2)
function _list_useless_plates(cuts, num_plates)
	reachable = falses(num_plates)
	reachable[1] = true
	num_iter_plate_removal = 0
	found_new_var = true
	while found_new_var
		found_new_var = false
		for (pp, fc, sc) in cuts
			if reachable[pp] && (!reachable[fc] || iszero(sc) || !reachable[sc])
				found_new_var = true
				reachable[fc] = true
				!iszero(sc) && (reachable[sc] = true)
			end
		end
		num_iter_plate_removal += 1
	end
	num_useless_plates = num_plates - sum(reachable)
	@show num_useless_plates
	@show num_iter_plate_removal
	return
end
=#

# TODO: check if this is used, and if it is, should be adapted to right
# notation.
function raw_warm_start(model, nzpe_idxs, nzpe_vals, nzcm_idxs, nzcm_vals)
	@assert length(nzpe_idxs) == length(nzpe_vals)
	@assert length(nzcm_idxs) == length(nzcm_vals)
	pe = model[:picuts]
	cm = model[:cuts_made]
	set_start_value.(pe[nzpe_idxs], nzpe_vals)
	set_start_value.(cm[nzcm_idxs], nzcm_vals)
end

# INTERNAL METHOD USED ONLY IN get_cut_pattern
# If the "pattern" is the extraction of a single piece return either:
# (1) a CutPattern representing the piece, if it is the same size as the
# original plate; otherwise (2) a CutPattern of the original plate containing
# a CutPattern representing the piece.
function extraction_pattern(
	bmr :: ByproductPPG2KP{D, S, P}, e_idx
) :: CutPattern{D, S} where {D, S, P}
	pli, pii = bmr.np[e_idx] # the plate index and the piece index
	L, W = bmr.pli_lwb[pli] # the plate dimensions
	piece = CutPattern(bmr.l[pii], bmr.w[pii], pii)
	L == bmr.l[pii] && W == bmr.w[pii] && return piece
	return CutPattern(L, W, false, CutPattern{D, S}[piece])
end

# INTERNAL METHOD USED ONLY IN get_cut_pattern
# Given a list of variables, returns two lists. The two lists have the same
# length (that is itself smaller than the given list). The first list is
# of indexes of the first list that have non-zero values (after aplying
# RoundNearest to them) and the second is the rounded (to nearest) values.
function gather_nonzero(vars, ::Type{D}) where {D}
	idxs = Vector{eltype(keys(vars))}()
	vals = Vector{D}()
	diff = 0.0
	for (idx, var) in pairs(vars)
		float_val = value(var)
		val = round(D, float_val, RoundNearest)
		abs(float_val - val) > diff && (diff = abs(float_val - val))
		iszero(val) && continue
		push!(idxs, idx)
		push!(vals, val)
	end
	diff > 1e-4 && @warn "At least one variable of a solved PPG2KP model had" *
		" a difference of $diff from the nearest integer. Maybe the tolerances" *
		" are too lax, or some variable that should not be relaxed is relaxed."
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
# SEE: https://discourse.julialang.org/t/unreachable-reached-at-0x7fa478093547-in-julia-1-0-5/36404
function build_cut_idx_stack(
	nz_cuts :: Vector{NTuple{3, P}},
	qt_cuts :: Vector{D},
	root_cut_idx :: Int #=P=#,
	debug :: Bool = false
) :: Vector{Int#=P=#} where {D, P}
	cut_idx_stack = Vector{Int#=P=#}()
	push!(cut_idx_stack, root_cut_idx)
	cuts_available = deepcopy(qt_cuts)
	@assert isone(cuts_available[root_cut_idx])
	cuts_available[root_cut_idx] -= one(D)
	next_cut = one(P)
	while next_cut <= length(cut_idx_stack)
		_, fc, sc = nz_cuts[cut_idx_stack[next_cut]]
		@assert !iszero(fc)
		fc_idx = consume_cut(nz_cuts, cuts_available, fc)
		fc_idx !== nothing && push!(cut_idx_stack, fc_idx)
		if !iszero(sc)
			sc_idx = consume_cut(nz_cuts, cuts_available, sc)
			sc_idx !== nothing && push!(cut_idx_stack, sc_idx)
		end
		next_cut += one(P)
	end
	# If cuts_available is different from a vector of zeros then some cuts
	# were not consumed (I am not sure if this should be possible).
	debug && @show cuts_available

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
	nzpe_idxs, nzpe_vals, bmr :: ByproductPPG2KP{D, S, P},
	debug :: Bool = false
) :: Dict{Int64, Vector{CutPattern{D, S}}} where {D, S, P}
	for (i, np_idx) in pairs(nzpe_idxs)
		pli, pii = bmr.np[np_idx]
		if debug
			@show np_idx
			@show pli, pii
			@show nzpe_vals[i]
			@show bmr.l[pii], bmr.w[pii]
		end
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
	nz_cut_idx_stack,
	nz_cuts :: Vector{NTuple{3, P}},
	nz_cuts_ori :: BitArray{1},
	bmr :: ByproductPPG2KP{D, S, P},
	debug :: Bool = false
) :: Dict{Int64, Vector{CutPattern{D, S}}} where {D, S, P}
	for cut_idx in reverse(nz_cut_idx_stack)
		debug && @show cut_idx
		pp, fc, sc = nz_cuts[cut_idx]
		debug && @show pp, fc, sc
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
		debug && @show ppl, ppw
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
	# local constant to alternate debug mode (method will not take a debug flag)
	debug = false
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

	# non-zero {piece extractions, cuts made} {indexes,values}
	nzpe_idxs, nzpe_vals = gather_nonzero(pe, D)
	nzcm_idxs, nzcm_vals = gather_nonzero(cm, D)
	if debug
		@show nzpe_idxs
		@show nzpe_vals
		@show nzcm_idxs
		@show nzcm_vals
		@show value.(pe[nzpe_idxs])
		@show value.(cm[nzcm_idxs])
	end

	sps_idx = check_if_single_piece_solution(bmr.np, nzpe_idxs)
	!iszero(sps_idx) && return extraction_pattern(bmr, sps_idx)

	# The cuts actually used in the solution.
	sel_cuts = bmr.cuts[nzcm_idxs]
	debug && @show sel_cuts
	# If the cut in `sel_cuts` is vertical or not.
	ori_cuts = nzcm_idxs .>= bmr.first_vertical_cut_idx
	# The index of the root cut (cut over the original plate) in `sel_cuts`.
	root_idx = findfirst(cut -> isone(cut[1]), sel_cuts)
	root_idx === nothing && return CutPattern(bmr.L, bmr.W, zero(D))

	cut_idx_stack = build_cut_idx_stack(sel_cuts, nzcm_vals, root_idx, debug)

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
	add_used_extractions!(patterns, nzpe_idxs, nzpe_vals, bmr, debug)

	bottom_up_tree_build!(patterns, cut_idx_stack, sel_cuts, ori_cuts, bmr, debug)

	if !isone(length(patterns)) || !haskey(patterns, 1) || !isone(length(patterns[1]))
		println("BUG AT GET_CUT_PATTERN")
		for (key, subpatts) in patterns
			@show key
			for subpatt in subpatts
				println(to_pretty_str(subpatt))
			end
		end
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
Note: `model[:plate_cons]` is a container of constraints of type and set
`(GenericAffExpr{Float64,VariableRef}, LessThan{Float64})` that represent all
the plate types (i.e., it will be the same length as ByproductPPG2KP.pli_lwb
and each position in `model[:plate_cons]` will correspond to the plate
described in the same position of `ByproductPPG2KP.pli_lwb`).
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
	end # @timeit "enumeration"

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

	@variable(model,
		0 <= cuts_made[i = 1:length(hvcuts)] <= pli_lwb[hvcuts[i][1]][3]
	, Int)

	# The objective function is to maximize the profit made by extracting
	# pieces from subplates.
	@objective(model, Max,
		sum(p[pii] * sum(picuts[pii2pair[pii]]) for pii = 1:num_piece_types)
	)

	# c1: There is just one of the original plate, and so it can be only used
	# to extract a single piece xor make a single cut that would make two new
	# subplates available.
	plate_number_one = @constraint(model,
		sum(picuts[pli2pair[1]]) + sum(cuts_made[parent2cut[1]]) <= 1
	)

	# c2: for each subplate type that is not the original plate, such subplate
	# type will be available the number of times it was the child of a cut,
	# subtracted the number of times it had a piece extracted or used for
	# further cutting.
	plate_numbers_2_to_m = @constraint(model, [ppli=1:(num_plate_types-1)],
		sum(picuts[pli2pair[ppli+1]]) + sum(cuts_made[parent2cut[ppli+1]]) <=
		sum(cuts_made[child2cut[ppli+1]])
	)
	# NOTE: we use a ppli that is pli-1 above because if we start the array at
	# 2 then we do not get an array but another container type. Also, these
	# constraints are not named, if we want to name them we can make a loop
	# below the next line, iterating the constraints and naming them.
	model[:plate_cons] = pushfirst!(plate_numbers_2_to_m, plate_number_one)

	if options["use-c25"]
		# c2.5: The amount of each subplate type generated by cuts (and used either
		# as a piece or as a intermediary plate) is bounded by the amount that can be
		# cut from the original plate.
		@constraint(model, c25[pli=2:num_plate_types],
			sum(picuts[pli2pair[pli]]) + sum(cuts_made[parent2cut[pli]]) <=
			pli_lwb[pli][3]
		)
	end

	# c3: the amount of each piece type extracted from different plate types
	# cannot surpass the demand for that piece type.
	@constraint(model,
		demand_con[pii=1:num_piece_types],
		sum(picuts[pii2pair[pii]]) <= d[pii]
	)

	lb = options["lower-bound"]
	if !iszero(lb)
		@constraint(model, obj_lb_con,
			sum(p[pii]*sum(picuts[pii2pair[pii]]) for pii = 1:num_piece_types) >= (lb + 1)
		)
	end

	ub = options["upper-bound"]
	if ub < sum(d .* p)
		@constraint(model, obj_ub_con,
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
	pp_dual = dual(constraints[pp])
	fc_dual = dual(constraints[fc])
	sc_dual = iszero(sc) ? 0.0 : dual(constraints[sc])
	#@show pp_dual
	#@show fc_dual
	#@show sc_dual
	# The reduced profit computation shown in the paper is done here.
	rp = fc_dual + sc_dual - pp_dual
	#@show rp
	return rp
end

# TODO: the call to the heuristic and solving the restricted MIP model
# probably should not be here, but on no_arg_check_build_model.
function _restricted_final_pricing!(
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
	#@timeit "relax_cm"
	cm_svcs = relax!(cm)
	#@timeit "relax_pe"
	pe_svcs = relax!(pe)
	#@assert length(cm) == length(bp.cuts)
	rc_vars = cm[rc_idxs]
	rc_svcs = cm_svcs[rc_idxs]
	# Every variable that is not a restricted cut is a unrestricted cut.
	uc_idxs = setdiff_sorted(keys(cm), rc_idxs)
	# Fix all unrestricted cuts (pool vars).
	uc_vars = cm[uc_idxs]
	#@timeit "fix_uc"
	fix.(uc_vars, 0.0; force = true)
	# Get the solution of the heuristic.
	bkv, _, shelves = Heuristic.iterated_greedy(
		d, p, bp.l, bp.w, bp.L, bp.W, rng
	)
	sol = Heuristic.shelves2cutpattern(shelves, bp.l, bp.w, bp.L, bp.W)
	LB = convert(Float64, bkv)
	# Solve the relaxed restricted model.
	flush_all_output()
	@timeit "lp_solve" optimize!(model)
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
	# TODO: if the assert below fails, then the solver used reaches an optimal
	# point of a LP but do not has duals for it, what should not be possible,
	# some problem has happened (for an example, the solver does not recognize
	# the model as a LP but think it is a MIP for example, CPLEX does this).
	@assert has_duals(model)
	#all_vars_values = value.(all_vars) # needs to be done before any changes
	rc_cuts = bp.cuts[rc_idxs]
	plate_cons = model[:plate_cons]
	debug && begin
		restricted_LB = LB
		restricted_LP = LP
		@show restricted_LB
		@show restricted_LP
	end

	# TODO: check if a tolerance is needed in the comparison. Query it from the
	# solver/model if possible (or use eps?).
	rc_discrete_bitstr = # the final pricing (get which variables should remain)
		@. floor(reduced_profit(rc_cuts, (plate_cons,)) + LP) >= LB
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
	#@timeit "fix_rc_subset"
	fix.(rc_vars[rc_fix_bitstr], 0.0; force = true)
	# The discrete variables are the restricted cuts that were not removed
	# by the final pricing of the restricted model and will be in the MIP
	# of the restricted model.
	discrete_vars = rc_vars[rc_discrete_bitstr]
	#@timeit "restore_ss_cm"
	restore!(discrete_vars, rc_svcs[rc_discrete_bitstr])
	# All piece extractions also need to be restored.
	#@timeit "restore_pe"
	restore!(pe, pe_svcs)
	#set_start_value.(all_vars, all_vars_values)
	flush_all_output()
	@timeit "mip_solve" optimize!(model) # restricted MIP solved
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
	#@timeit "relax_ss_rc"
	relax!(rc_vars[rc_discrete_bitstr])
	#@timeit "relax_pe"
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
function _recompute_idxs_to_add!(
	pool_idxs_to_add, pool, plate_cons, threshold, n_max :: P,
) :: Tuple{P, Bool} where {P}
	found_above_threshold = false
	num_positive_rp_vars = zero(P)
	# TODO: consider if the optimization of using the pool_idxs_to_add as a
	# buffer and deciding between overwrite and push! is worth.
	empty!(pool_idxs_to_add)

	for (pool_idx, cut) in pairs(pool)
		rp = reduced_profit(cut, plate_cons)
		rp > 0.0 && continue # non-positive reduced profit is irrelevant to us
		num_positive_rp_vars += one(P)
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

	# Only n_max variables are added/unfixed in a single iteration.
	n_max <= length(pool_idxs_to_add) && resize!(pool_idxs_to_add, n_max)

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
function _iterative_pricing!(
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
	#=
	if JuMP.solver_name(model) == "Gurobi"
		old_method = get_optimizer_attribute(model, "Method")
		# Change the LP-solving method to dual simplex, this allows for better
		# reuse of the partial solution.
		set_optimizer_attribute(model, "Method", 1)
		set_optimizer_attribute(model, "Presolve", 0)
	end
	=#
	flush_all_output()
	# the last solve before this was MIP and has no duals
	@timeit "solve_lp" optimize!(model)
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
	resize!(to_unfix, n_max)
	# Not sure if restarting the model use the old values as start values.
	#all_vars = all_variables(model)
	num_positive_rp_vars = initial_num_positive_rp_vars
	# The iterative pricing continue until there are variables to unfix.
	while !isempty(to_unfix)
		debug && begin
			@show num_positive_rp_vars
			@show length(to_unfix)
			@show was_above_threshold
		end
		@timeit "update_bounds" begin
			unfixed_vars = unused_vars[to_unfix]
			unfix.(unfixed_vars)
			set_lower_bound.(unfixed_vars, unused_lbs[to_unfix])
			set_upper_bound.(unfixed_vars, unused_ubs[to_unfix])
		end
		@timeit "deleteat!" begin
			deleteat!(unused_vars, to_unfix)
			deleteat!(unused_cuts, to_unfix)
			deleteat!(unused_lbs, to_unfix)
			deleteat!(unused_ubs, to_unfix)
		end
		@assert allsame(length.((unused_vars, unused_cuts, unused_lbs, unused_ubs)))
		#set_start_value.(all_vars, value.(all_vars))
		flush_all_output()
		@timeit "solve_lp" optimize!(model)
		flush_all_output()
		@timeit "recompute" begin
			num_positive_rp_vars, was_above_threshold = _recompute_idxs_to_add!(
				to_unfix, unused_cuts, plate_cons, threshold, n_max
			)
		end
		was_above_threshold && (pricing_threshold_hits += 1)
	end

	debug && @show pricing_threshold_hits
	#=
	if JuMP.solver_name(model) == "Gurobi"
		# Undo the change done before (look at where old_method is defined).
		set_optimizer_attribute(model, "Method", old_method)
		set_optimizer_attribute(model, "Presolve", -1)
	end
	=#

	return
end

function _partition_by_bits(bits, list)
    accepted = Vector{eltype(list)}(undef, length(list))
    rejected = Vector{eltype(list)}(undef, length(list))
    qt_accepted = qt_rejected = 0
    for (istrue, value) in zip(bits, list)
        if istrue
            qt_accepted += 1
            accepted[qt_accepted] = value
        else
            qt_rejected += 1
            rejected[qt_rejected] = value
        end
    end
    resize!(accepted, qt_accepted)
    resize!(rejected, qt_rejected)
    return accepted, rejected
end

# TODO: check if the kept variables are chosen from the unfixed variables
# of from all variables.
function _final_pricing!(
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
	#println(string(floor.(reduced_profit.(bp.cuts, (plate_cons,)))))
	to_keep_bits = # the final pricing (get which variables should remain)
		@. floor(reduced_profit(bp.cuts, (plate_cons,)) + LP) >= LB
	#to_keep_idxs = key(vars)[unfixed_bits][to_keep_bits]
	to_keep_idxs, to_delete_idxs = _partition_by_bits(to_keep_bits, keys(vars))
	new_fvci = searchsortedfirst(
		to_keep_idxs, bp.first_vertical_cut_idx
	)
	# The two @show below are too verbose to be useful except when debugging
	# small instances.
	#@show bp.cuts[to_delete_idxs]
	#@show bp.cuts[to_keep_idxs]
	debug && begin
		final_pricing_num_vars_removed = length(to_delete_idxs)
		final_pricing_num_vars_kept = length(to_keep_idxs)
		@show final_pricing_num_vars_removed
		@show final_pricing_num_vars_kept
	end
	deleteat!(bp.cuts, to_delete_idxs)
	to_keep_vars, to_delete_vars = _partition_by_bits(to_keep_bits, vars)
	JuMP.delete(model, to_delete_vars)
	model[:cuts_made] = to_keep_vars
	new_bp = ByproductPPG2KP(
		bp.cuts, new_fvci, bp.np, bp.pli_lwb, bp.l, bp.w, bp.L, bp.W
	)
	return new_bp, to_keep_bits
end

function _furini_pricing!(
	model, byproduct, d, p, seed, alpha, beta, debug
) where {D, S, P}
	rng = Xoroshiro128Plus(seed)
	rc_idxs = all_restricted_cuts_idxs(byproduct)
	@timeit "restricted" bkv, sol, pe_svcs, cm_svcs = _restricted_final_pricing!(
		model, rc_idxs, d, p, byproduct, rng, debug
	)
	UB = convert(Float64, bkv)
	# If those fail, we need may need to rethink the iterative pricing,
	# because the use of max_profit seems to assume no negative profit items.
	@assert all(e -> e >= zero(e), d)
	@assert all(e -> e >= zero(e), p)
	# TODO: check if is needed to pass the whole bp structure
	@timeit "iterative" _iterative_pricing!(
		model, byproduct.cuts, cm_svcs, sum(d .* p), alpha, beta, debug
	)
	LB = objective_value(model)
	# TODO: change the name of the methods to include the exclamation mark
	@timeit "final" byproduct, to_keep_bits = _final_pricing!(
		model, byproduct, UB, LB, debug
	)
	debug && @timeit "plate_check" begin
		total_num_plates = length(byproduct.pli_lwb)
		num_r, r = _reachable_plate_types(byproduct.cuts, total_num_plates)
		@assert num_r <= total_num_plates
		@assert length(r) == total_num_plates
		num_unreachable_plates = total_num_plates - num_r
		@show num_unreachable_plates
		println("START_UNREACHABLE_PLATES")
		if num_unreachable_plates > 0
			for (plate_idx, isreachable) in zip(1:total_num_plates, r)
				if !isreachable
					l, w, b = byproduct.pli_lwb[plate_idx]
					println("i $plate_idx l $l w $w b $b")
				end
			end
		end
		println("END_UNREACHABLE_PLATES")
	end

	# Restore the piece extractions and the kept variables to their original
	# configuration. Note that cuts_made is now a subset of what it was before.
	@timeit "last_restore" @views restore!(
		[model[:picuts]; model[:cuts_made]], [pe_svcs; cm_svcs[to_keep_bits]]
	)

	return byproduct
end

# SEE SECTION 4.2: unfortunately they decided to complicate the method
# even more adding two parameters alpha and beta.
# List items 1 and 3 (not numbered) of the first list of section 4.4
# complicate things further. In fact, both the greedy heuristic and the
# restricted model are solved but with a time limit (the best solution
# found in the middle of this process is used for the true final pricing).
function _no_arg_check_build_model(
	model, d :: Vector{D}, p :: Vector{P}, l :: Vector{S}, w :: Vector{S},
	L :: S, W :: S, options :: Dict{String, Any} = Dict{String, Any}()
) :: ByproductPPG2KP{D, S, P} where {D, S, P}
	byproduct = build_complete_model(model, d, p, l, w, L, W, options)
	debug = options["verbose"] && !options["quiet"]
	debug && begin
		length_pe_before_pricing = length(model[:picuts])
		length_cm_before_pricing = length(model[:cuts_made])
		@show length_pe_before_pricing
		@show length_cm_before_pricing
	end
	pricing_method = options["pricing"]
	pricing_method == "none" && return byproduct

	@timeit "pricing" begin
		if pricing_method == "expected"
			pricing_method = (options["faithful2furini2016"] ? "furini" : "becker")
		end

		if pricing_method == "furini"
			_furini_pricing!(
				model, byproduct, d, p, options["pricing-heuristic-seed"],
				options["pricing-alpha"], options["pricing-beta"], debug
			)
		else
			@assert pricing_method == "becker"
			@warn "The 'becker' pricing was not yet implemented. Acting as 'none' was passed."
		end
	end # @timeit "pricing"

	return byproduct
end

# Only used to define `build_model(Val{:PPG2KP}, ...)`.
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
	return _no_arg_check_build_model(model, d, p, l, w, L, W, norm_options)
end

end # module

