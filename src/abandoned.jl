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
