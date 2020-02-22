module PPG2KP
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

# TODO: document
# TODO: think of better names for the fields? maybe update the names on the
#       model-building code?
# TODO: check the abandoned.jl to see if parts of the original instance are
#       needed to create the solution
# TODO: maybe add this to Enumeration, use as return of gen_cuts and re-export
#       in this module
struct ByproductPPG2KP{D, S, P}
	hcuts :: Vector{NTuple{3, P}}
	vcuts :: Vector{NTuple{3, P}}
	np :: Vector{Tuple{P, D}}
	pli_lwb :: Vector{Tuple{S, S, P}}
	faithful2furini2016 :: Bool
end

include("./Heuristic.jl")
include("./Args.jl")

include("./Enumeration.jl")
using .Enumeration
using ..Utilities

using JuMP
using TimerOutputs

import ..get_cut_pattern, ..CutPattern
function get_cut_pattern(
	model_type :: Val{:PPG2KP}, model :: JuMP.Model, ::Type{D}, ::Type{S},
	build_model_return :: ByproductPPG2KP{D, S}
) :: CutPattern{D, S} where {D, S}
	bmr = build_model_return
	# 1. Check if there can be an extraction from the original plate to a
	#    single piece. If it may and it happens, then just return this single
	#    piece solution; otherwise the first cut is in `cuts_made`, find it
	#    (i.e., traverse non-zero cuts_made variables checking for a hcut or
	#    vcut that has the original plate as the parent plate).
	# 2. Given the children of the first cut, find them in either `picuts`
	#    or `cuts_made`, if they are in `picuts` they are an extraction (
	#    this may be a plate into a piece, or immediately a piece, in the
	#    case the plate has the exact size of the piece); if they are in
	#    `cuts_made` we need to discover the orientation (find in {h,v}cuts)
	#    create a CutPattern and add to a list of pending branches to open;
	#    if they are not in either `picuts` or `cuts_made` then they
	#    are waste and a leaf plate may be added.
	#=pe_ = model[:picuts] # Piece Extractions
	# USE: https://docs.julialang.org/en/v1.0/stdlib/SparseArrays/#Sparse-Vector-Storage-1
	cm = model[:cuts_made]
	ps_nz_idx = (i in keys(ps) if is_valid(m, ps[i]) && value(ps[i]) ≉ 0.0)
	cm_nz_idx = (i in keys(cm) if is_valid(m, cm[i]) && value(cm[i]) ≉ 0.0)

	for to_
	for @view vcat(bmr.vnnn, bmr.hnnn)
	end
	println("(plate length, plate width) => (number of times this extraction happened, piece length, piece width)")
	foreach(ps_nz) do e
		i, v = e
		pli, pii = np[i]
		println((pli_lwb[pli][1], pli_lwb[pli][2]) => (v, l[pii], w[pii]))
	end
	println("(parent plate length, parent plate width) => (number of times this cut happened, (first child plate length, first child plate width), (second child plate length, second child plate width))")
	foreach(cm_nz) do e
		i, v = e
		parent, fchild, schild = hvcuts[i]
		@assert !iszero(parent)
		@assert !iszero(fchild)
		if iszero(schild)
			@assert p_args["faithful2furini2016"]
			if pli_lwb[parent][1] == pli_lwb[fchild][1]
							println((pli_lwb[parent][1], pli_lwb[parent][2]) => (v, (pli_lwb[fchild][1], pli_lwb[fchild][2]), (pli_lwb[parent][1], pli_lwb[parent][2] - pli_lwb[fchild][2])))
			else
				@assert pli_lwb[parent][2] == pli_lwb[fchild][2]
							println((pli_lwb[parent][1], pli_lwb[parent][2]) => (v, (pli_lwb[fchild][1], pli_lwb[fchild][2]), (pli_lwb[parent][1] - pli_lwb[fchild][1], pli_lwb[parent][2])))
			end
		else
						println((pli_lwb[parent][1], pli_lwb[parent][2]) => (v, (pli_lwb[fchild][1], pli_lwb[fchild][2]), (pli_lwb[schild][1], pli_lwb[schild][2])))
		end
		end
	end=#
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
    disable_unrestricted_cuts!(m, sl, sw, nnn, pli_lwb)

!!! **Internal use.**

"""
function disable_unrestricted_cuts!(m, sl, sw, nnn, pli_lwb)
	@assert length(sl) == length(sw)
	n = length(sl)
	@assert issorted(sl)
	@assert issorted(sw)
	reg = Vector{SavedBound}()
	# For each triple in nnn, there is an associated variable in the
	# cuts_made vector of m. The orientation and the position of the cut
	# may be obtained by: getting the parent plate and first child indexes
	# from nnn, using them to index pli_lwb, check which dimension has
	# been changed and its new size. If a vertical (horizontal) cut creates
	# a first child with size s, and s is NOT present in sw (sl), then that
	# varible will be disabled.
	for (i ,(pp, fc, _)) in enumerate(nnn)
		ppl, ppw, _ = pli_lwb[pp]
		fcl, fcw, _ = pli_lwb[fc]
		@assert fcl < ppl || fcw < ppw
		should_fix = false
		if fcl < ppl
			fcl_idx = searchsortedfirst(sl, fcl)
			if fcl_idx > n || sl[fcl_idx] != fcl
				should_fix = true
			end
		else
			@assert fcw < ppw
			fcw_idx = searchsortedfirst(sw, fcw)
			if fcw_idx > n || sw[fcw_idx] != fcw
				should_fix = true
			end
		end
		if should_fix
			var = m[:cuts_made][i]
			# fix(...; force = true) erase the old bounds to then fix the variable.
			# To be able to restore variables that had bounds, it is necessary to
			# save the old bounds and the varible reference to restore them.
			save_bound_if_exists!(reg, var)
			fix(var, 0.0; force = true)
		end
	end
	reg
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
function no_arg_check_build_model(
	model, d :: Vector{D}, p :: Vector{P}, l :: Vector{S}, w :: Vector{S},
	L :: S, W :: S, options :: Dict{String, Any} = Dict{String, Any}()
) where {D, S, P}
	@timeit "enumeration" begin
	@assert length(d) == length(l) && length(l) == length(w)
	num_piece_types = convert(D, length(d))

	sllw = SortedLinkedLW(D, l, w)
	# TODO: change enumeration to also use a Dict?
	pli_lwb, hcuts, vcuts, np = gen_cuts(P, d, sllw, L, W;
		ignore_2th_dim = options["ignore-2th-dim"],
		ignore_d = options["ignore-d"],
		round2disc = options["round2disc"],
		no_cut_position = options["no-cut-position"],
		no_redundant_cut = options["no-redundant-cut"],
		no_furini_symmbreak = options["no-furini-symmbreak"],
		faithful2furini2016 = options["faithful2furini2016"]
	)
	num_plate_types = length(pli_lwb)
	hvcuts = vcat(hcuts, vcuts)

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
		@variable(model, picuts[1:length(np)] >= 0, Int)
		#@variable(model,
		#  0 <= picuts[i = 1:length(np)] <=
		#    min(pli_lwb[np[i][1]][3], d[np[i][2]]),
		#Int)
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

	model, hvcuts, pli_lwb, np
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
) where {D, S, P}
	norm_options = create_normalized_arg_subset(
		options, accepted_arg_list(Val(:PPG2KP))
	)
	no_arg_check_build_model(model, d, p, l, w, L, W, norm_options)
	return model
end

end # module

