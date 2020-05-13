# This file is not included by any module and is not intended to exist in
# the long-term. It is just a home for julia code that was developed and
# needs to be renamed, moved to a new module, or maybe deleted.

function build(
	model, ::Type{P}, d :: Vector{D}, l :: Vector{S}, w :: Vector{S},
	L :: S, W :: S
) where {D, S, P}
	@assert length(l) == length(w) && length(w) == length(d)
	a = l .* w # area of the piece types

	partitions(P, d, l, w, L, W)
	@variables model begin
		x[1:N, 1:T], Bin  # true if piece type is placed at node, false otherwise
		v[1:N, 1:DL], Bin # true if node is a vertical cut, false otherwise
		h[1:N, 1:DW], Bin # true if node is a horizontal cut, false otherwise
		l_[1:N] >= 0 # length available to node N
		w_[1:N] >= 0 # width available to node N
	end

	@objective(model, Max, sum(p[j] * x[i, j] for i = 1:N, j = 1:T))

	# ub0: trivial area upper bound
	@constraint(model, sum(a[j] * x[i, j] for i = 1:N, j = 1:T) <= L*W)

	# c1: Each node is either: 1) a placed piece (leaf); 2) a horizontal cut
	# (branch); 3) a vertical cut (branch).
	for i = 1:N
		@constraint(model,
			sum(x[i, j] for j = 1:T) + sum(v[i, j] for j = 1:DL) +
			sum(h[i, j] for j = 1:DW) <= 1
		)
	end
	# c2: Each non-root node can only be a branch or leaf if the node parent
	# is a branch.
	for i = 2:N
		@constraint(model,
			sum(x[i, j] for j = 1:T) + sum(v[i, j] for j = 1:DL) +
			sum(h[i, j] for j = 1:DW) <= sum(v[div(i, 2), j] for j = 1:DL) +
			sum(h[div(i, 2), j] for j = 1:DW)
		)
	end

	# c3 and c4: The root node is limited by the size of the original plate.
	fix(l_[1], L; force = true)
	fix(w_[1], W; force = true)

	# c5: There is a limit on the number of pieces available for each piece type.
	for j = 1:T
		@constraint(model,
			sum(x[i, j] for i = 1:N) <= d[j]
		)
	end

	# c6 and c7: The cuts reduce the space available for their children.
	# c9 and c10: the length (width) of a plate is always at max the length
	# NOTE: it is not possible to tighten the big-M using dl[1] and dw[1] because
	# they can be children that yet keep one dimension at the same size as root.
	for i = 2:N
		i2 = div(i, 2)
		if iseven(i) # the first child needs a big-M to restrict its dimensions
			@constraints model begin
				l_[i] <= sum(dl[j] * v[i2, j] for j = 1:DL) + L * sum(h[i2, j] for j = 1:DW)
				w_[i] <= sum(dw[j] * h[i2, j] for j = 1:DW) + W * sum(v[i2, j] for j = 1:DL)
				# these are c9 and c10, they are implicit for the isodd(i) case
				l_[i] <= l_[i2]
				w_[i] <= w_[i2]
			end
		else
			# the second child does not need a big-M to restrict its dimensions,
			# nor does it need c9 and c10
			@constraints model begin
				l_[i] <= l_[i2] - sum(dl[j] * v[i2, j] for j = 1:DL)
				w_[i] <= w_[i2] - sum(dw[j] * h[i2, j] for j = 1:DW)
			end
		end
	end
	# c11 and c12: the placed pieces respect length (width) of the plate where
	# they are placed.
	for i = 1:N
		@constraints model begin
			sum(l[j] * x[i, j] for j = 1:T) <= l_[i]
			sum(w[j] * x[i, j] for j = 1:T) <= w_[i]
		end
	end

	model
end # build

# HIGH LEVEL EXPLANATION OF THE MODEL
#
# Variables:
#
# `pieces_sold[pii]`: Integer. The number of pieces `pii` sold. Is the minimum
#   between the demand of the piece and the amount of plates generated and not
#   used that have exactly the same size as the piece.
# `cuts_made[n1, n2, n3]`: Integer. The number of subplates of type `n1` that
#   are cut into subplates `n2` and `n3` (horizontal and vertical cuts are
#   together for now). As this is a symmetry-breaking model, the plate types
#   are not only each distinct `l` and `w` but each different `l`, `w`, and
#   `symm` (that marks if the plate can be cut only horizontally, only
#   vertically, or both ways).
#
# Objective function:
#
# Maximize the profit of the pieces sold.
#   sum(p[pii] * pieces_sold[pii])
#
# Constraints:
#
# There is exactly one of the original plate, which may be used for cutting
# or extracting a piece.
#   sum(pieces_sold[plates exactly the size of the original plate]) +
#   sum(cuts_made[1, _,  _]) <= 1
# The number of subplates available depends on the number of plates that have
# it as children.
#   sum(cuts_made[n1>1, _, _]) <= sum(cuts_made[_, n2, n3])
#     where n2 == n1 or n3 == n1, doubling cuts_made[_, n2, n3] if n2 == n3
# The number of pieces sold is bounded both by the demand of the piece type and
# the the number of unused plates with the same size as the piece.
#   sum(pieces_sold[pii]) <= d[pii]
#   sum(pieces_sold[pii]) <= cuts_made[_, n2, n3] - cuts_made[n2 or n3, _, _]
#     where n2 or n3 has the same size as pii, fixing for when n2 == n3
#
# Unnecessary constraints:
#
# The amount of times a plate type may be cut is bounded by how many of them
# could fit the original plate. Note that we ignore the symmetry tag here
# and group all the plates with the same `l` and `w` but distinct symmetry tag.
#   sum(cuts_made[plates sharing `l` and `w`, _, _]) <= (L ÷ l) * (W ÷ w)
function build_model_with_symmbreak(
	model, d :: Vector{D}, p :: Vector{P}, l :: Vector{S}, w :: Vector{S},
	L :: S, W :: S; only_binary = false, use_c25 = false,
	ignore_2th_dim = false, ignore_d = false, round2disc = true
) where {D, S, P}
	@assert length(d) == length(l) && length(l) == length(w)
	num_piece_types = convert(D, length(d))

	sllw = SortedLinkedLW(D, l, w)
	pli2lwsb, hcuts, vcuts, pii2plis, pli2piis, same_size_plis =
		gen_cuts_sb(P, d, sllw, L, W; ignore_2th_dim = ignore_2th_dim,
		ignore_d = ignore_d,
		round2disc = round2disc
	)
	num_plate_types = length(pli2lwsb)
	hvcuts = vcat(hcuts, vcuts)

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
		!iszero(schild) && push!(child2cut[schild], i)
	end

	# If all pieces have demand one, a binary variable will suffice to make the
	# connection between a piece type and the plate it is extracted from.
	naturally_only_binary = all(di -> di <= 1, d)
	if naturally_only_binary || only_binary
		# only_binary is equal to naturally_only_binary for now, but the idea is
		# that only_binary will expand the number of binary variables to account
		# for piece solds that can repeat (i.e., have demand more than one).
		# NOTE that using only_binary with this model will restrict much more
		# than using only_binary with the model without symmetry, unless the demand
		# of all pieces is naturally_only_binary the model will give much worse
		# results.
		@variable(model, pieces_sold[1:num_piece_types], Bin)
	else
		@variable(model,
			0 <= pieces_sold[pii = 1:num_piece_types] <= d[pii],
		Int)
	end

	if only_binary
		@variable(model, cuts_made[1:length(hvcuts)], Bin)
	else
		@variable(model, cuts_made[1:length(hvcuts)] >= 0, Int)
	end

	# The objective function maximizes the profit of the pieces sold.
	@objective(model, Max,
		sum(p[pii] * sum(pieces_sold[pii]) for pii = 1:num_piece_types)
	)

	# c1: There is just one of the original plate, and so it can be only used to
	# extract a single piece (this is rare, because there would need to be a
	# piece the exact same size as the original plate) xor make a single cut that
	# would make two new subplates available.
	@constraint(model,
		sum(pieces_sold[pli2piis[1]]) + sum(cuts_made[parent2cut[1]]) <= 1
	)

	# c2: for each subplate type that is not the original plate, such subplate
	# type may be cut at most the number of times it was a child of another cut.
	for pli in 2:num_plate_types
		@constraints model begin
			sum(cuts_made[parent2cut[pli]]) <= sum(cuts_made[child2cut[pli]])
		end
	end

	if use_c25
		@error "not sure if correctly implemented, check before"
		# TODO: check if the below is corret. What seems to be wrong is that ssplis
		# is a vector of plate indexes, while cuts_made should be indexed by
		# vectors of cut indexes.

		# c2.5: The amount of each subplate type generated by cuts (and used either
		# as a piece or as a intermediary plate) is bounded by the amount that can be
		# cut from the original plate.
		for ssplis in same_size_plis
			@assert !iszero(length(ssplis))
			@assert isone(length(unique!(map(i -> pli2lwsb[i][1], ssplis))))
			@assert isone(length(unique!(map(i -> pli2lwsb[i][2], ssplis))))
			@constraint(model,
				sum(cuts_made[ssplis]) <= pli2lwsb[ssplis[1]][4]
			)
		end
	end

	# c3: finally, for each piece type, the amount of pieces sold of that type is
	# at most the number of plates with the piece exact size that were not cut to
	# make smaller plates.
	for pii in 1:num_piece_types
		@constraint(model, pieces_sold[pii] <= sum(cuts_made[vcat(child2cut[pii2plis[pii]]...)]) - sum(cuts_made[vcat(parent2cut[pii2plis[pii]]...)]))
	end

	model, hvcuts, pli2lwsb, pii2plis, pli2piis
end # build_model_with_symmbreak


#=
# CODE that lived inside the script method, and printed the solutions
# from the models. Need to be broken in many model-specific methods
# that take a solved model, and maybe the model building output, and
# create a common and structured representation of the solution.
		if p_args["flow-model"]
			ps = value.(m[:edge][pii] for pii = 1:length(d) if is_valid(m, m[:edge][pii]))
			ps_nz = [iv for iv in enumerate(ps) if iv[2] > 0.001]
			@show ps_nz
		else
			if p_args["break-hvcut-symmetry"]
				ps = m[:pieces_sold]
				cm = m[:cuts_made]
				ps_nz = [(i, value(ps[i])) for i = 1:length(ps) if is_valid(m, ps[i]) && value(ps[i]) > 0.001]
				cm_nz = [(i, value(cm[i])) for i = 1:length(cm) if is_valid(m, cm[i]) && value(cm[i]) > 0.001]
				@show ps_nz
				@show cm_nz
				println("(piece length, piece width) => (piece index, amount in solution, profit of single piece, total profit contributed in solution) ")
				foreach(ps_nz) do e
					pii, v = e
					println((l[pii], w[pii]) => (pii, v, p[pii], v * p[pii]))
				end
				println("(parent plate length, parent plate width) => ((first child plate length, first child plate width), (second child plate length, second child plate width))")
				foreach(cm_nz) do e
					i, _ = e
					parent, fchild, schild = hvcuts[i]
					if iszero(schild)
						println((pli2lwsb[parent][1], pli2lwsb[parent][2]) => ((pli2lwsb[fchild][1], pli2lwsb[fchild][2]), (0, 0)))
					else
						println((pli2lwsb[parent][1], pli2lwsb[parent][2]) => ((pli2lwsb[fchild][1], pli2lwsb[fchild][2]), (pli2lwsb[schild][1], pli2lwsb[schild][2])))
					end
				end
			else # same as: if !p_args["break-hvcut-symmetry"]
				p_args["relax2lp"] && @warn "relax2lp flag used, be careful regarding solution"
				ps_nz = Vector{Tuple{Int, Float64}}()
				cm_nz = Vector{Tuple{Int, Float64}}()
				for i = 1:length(ps)
					if is_valid(m, ps[i]) && value(ps[i]) > 0.0001
						push!(ps_nz, (i, value(ps[i])))
					end
				end
				for i = 1:length(cm)
					if is_valid(m, cm[i]) && value(cm[i]) > 0.0001
						push!(cm_nz, (i, value(cm[i])))
					end
				end
				#ps_nz = [(i, value(ps[i])) for i = 1:length(ps) if (is_valid(m, ps[i]) && value(ps[i]) > 0.001)]
				#cm_nz = [(i, value(cm[i])) for i = 1:length(cm) if (is_valid(m, cm[i]) && value(cm[i]) > 0.001)]
				@show ps_nz
				@show cm_nz
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
			end
		end
	end
=#

# A situation that occurred with some LP solving methods was the following:
# The method below intends to fix an unproductive behavior that happened with
# some LP methods. The method would give a very small positive reduced profit
# for a few variables, and the _iterative_pricing would keep executing because
# of it. The values were so small and the variables so few that they could not
# affect the objective value of the LP significantly, but triggered a
# resolving anyway.
function _is_significant(LP, cuts, plate_cons)
	delta_to_next_int = (one(LP) + floor(LP)) - LP
	improvement = sum(_reduced_profit.(cuts, (plate_cons,)))
	@show improvement
	@show delta_to_next_int
	return improvement >= delta_to_next_int
end


#################### START OF THE OLD WARM_START FILE #################### 
"""
    _search_approx_cut(pp, fcl, fcw, max_diff, approx_l, nnn, pli_lwb) :: P

!!! **Internal use.**

"""
function _search_approx_cut(
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
    _search_cut_or_symmetry(pp, fcl, fcw, nnn, pli_lwb)

!!! **Internal use.**

"""
function _search_cut_or_symmetry(
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
    _search_cut(pp, fcl, fcw, nnn, pli_lwb)

!!! **Internal use.**

"""
function _search_cut(
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

# TODO: Check why this method is used if a structure like SortedLinkedLW
# would answer this more efficiently and be aware of the type used for the
# index.
"""
    _min_l_fitting_piece(l, w, L, W)

!!! **Internal use.**

Given a plate `L`x`W` and two pieces dimensions list (`l`, `w`),
return the index of the piece of smallest length that fits the
plate (the width dimension may preclude this from being just
the piece of smallest length).
"""
function _min_l_fitting_piece(l, w, L, W)
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

# Warm start faithful2furini.
function warm_start_f2f(
	model, l, w, L, W,
	pat :: Vector{Vector{D}},
	bp :: ByproductPPG2KP{D, S, P}
	#round2disc wait to see if this is needed
	# which other model building options will need to be passed to this?
) where {D, S, P}
	# TODO: implement
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
			cut = _search_cut(rpli, L, ws, nnn, pli_lwb)
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
				cut = _search_cut(rpli, lp, ws, nnn, pli_lwb)
				cut_var_vals[cut] = 1 + get!(cut_var_vals, cut, 0)
				_, fc, rpli = nnn[cut]

				fcl, fcw, _ = pli_lwb[fc]
				println("$(i)\ti\t$(fcl)\t$(fcw)")
				println("$(piece)\tp\t$(lp)\t$(wp)")

				min_w = _min_l_fitting_piece(w, l, fcw, fcl)
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
						cut2 = _search_cut(fc, lp, wp, nnn, pli_lwb)
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
						cut2 = _search_approx_cut(
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
				min_l_pii = _min_l_fitting_piece(l, w, pll, plw)
				min_l = iszero(min_l_pii) ? pll : l[min_l_pii]
				if pll >= lp + min_l
					triml = pll - lp
					cut = _search_approx_cut(
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
				min_w_pii = _min_l_fitting_piece(w, l, plw, pll)
				min_w = iszero(min_w_pii) ? plw : l[min_w_pii]
				if plw >= wp + min_w
					trimw = plw - wp
					cut = _search_approx_cut(
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

# TODO: Consider if this should be enabled back again considering that, for
# now, the idea of a symmetry-breaking version of PP-G2KP was abandoned.
# TODO: internal change. The symmetry tag dimension of plis should be the
# first one, given how columnwise memory layout works.
#=
function gen_cuts_sb(
	::Type{P}, d :: Vector{D}, sllw :: SortedLinkedLW{D, S}, L :: S, W :: S;
	ignore_2th_dim = false, ignore_d = false, round2disc = true
) where {D, S, P}
	(ignore_d || ignore_2th_dim || round2disc) && @error "ignore_2th_dimm, ignore_d, and round2disc are not yet implemented for gen_cuts_sb, also, first improve the performance by memoizing the discretizations"
	l = sllw.l
	w = sllw.w
	@assert length(d) == length(l)
	@assert length(d) == length(w)
	only_single_l = discretize(
		d, l, w, L, W; only_single_pieces = true, ignore_W = ignore_2th_dim,
		ignore_d = ignore_d
	)
	only_single_w = discretize(
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
			for y in discretize(
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
			for x in discretize(
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
=#

