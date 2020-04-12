module PPG2KP
# 2.2) Finally, check which variables are fixed to zero and remove them both
#   from the model as from the ByproductPPG2KP. We need to check if we will
#   remove constraints that become irrelevant, if this is done we will need
#   to update all the cuts to refer to the new plate indexes.
# Include submodules.
include("./Heuristic.jl")
include("./Args.jl")
include("./Enumeration.jl")

using .Enumeration
export ByproductPPG2KP # re-export ByproductPPG2KP from Enumeration

using ..Utilities
using RandomNumbers.Xorshifts # for Xoroshiro128Plus (for heuristic rng)
import ..to_pretty_str # for debug

# Import third-party libraries.
using JuMP
using TimerOutputs
import MathOptInterface
const MOI = MathOptInterface

# Include code that is focused on implementing a specific feature, but does
# not merit a module of their own. Non-exported methods of these files are not
# expected to be used in this module, but exported methods are exports from
# this module (ou specialization of parent modules).
include("./get_cut_pattern.jl")
include("./warm_start.jl")

function _partition_by_bits(bits, list)
	@assert length(bits) == length(list) # yes, it is length, not axes, what we want
	accepted = similar(list)
	rejected = similar(list)
	indexes = LinearIndices(list)
	qt_accepted = qt_rejected = 0 # yes, this does not need to be a specific type
	for (istrue, value) in zip(bits, list)
		if istrue
			qt_accepted += 1
			accepted[indexes[qt_accepted]] = value
		else
			qt_rejected += 1
			rejected[indexes[qt_rejected]] = value
		end
	end
	resize!(accepted, qt_accepted)
	resize!(rejected, qt_rejected)
	return accepted, rejected
end

# Delete the vars from the model and from the byproduct.
# Returns byproduct because if it is :cuts_made, then the old object is
# invalid and a new one is necessary (the vectors are changed inplace, but
# the first_vertical_cut_idx cannot be updated because the struct is
# immutable). If it is :picuts, then the passed struct is updated
# to a valid state and then returned.
function _delete_vars!(bp, model, to_keep_bits, var_set_name_in_model)
	@assert var_set_name_in_model in (:picuts, :cuts_made)
	var_set_name_in_bp = var_set_name_in_model == :picuts ? :np : :cuts
	vars = model[var_set_name_in_model]
	to_keep_idxs, to_delete_idxs = _partition_by_bits(to_keep_bits, keys(vars))
	deleteat!(getfield(bp, var_set_name_in_bp), to_delete_idxs)
	to_keep_vars, to_delete_vars = _partition_by_bits(to_keep_bits, vars)
	JuMP.delete(model, to_delete_vars)
	model[var_set_name_in_model] = to_keep_vars
	var_set_name_in_model == :picuts && return bp
	new_fvci = searchsortedfirst(to_keep_idxs, bp.first_vertical_cut_idx)
	return ByproductPPG2KP( # new_fvci was computed and need to be updated
		bp.cuts, new_fvci, bp.np, bp.pli_lwb, bp.l, bp.w, bp.L, bp.W
	)
end

# Include all code related to finding and removing unreachable plates/vars.
# It does make use of _partition_by_bits and _delete_vars! but both are
# used elsewhere so they are kept in this file. From the methods defined
# in the include below only _purge_unreachable! is used here.
include("./unreachable.jl")

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

	# TODO: Optional. Consider if changing the order of the variables can impact
	# the solver performance. The only problem with this is that current
	# ByproductPPG2KP forces all horizontal cuts to come before all vertical
	# cuts. Changing the order of the constraints, however, needs the complete
	# recomputation of ByproductPPG2KP.cuts and ByproductPPG2KP.np.
	# Some criteria to be considered for the np variables (mostly):
	# how much is the absolute profit obtained; how much relative area is wasted;
	# how much is the profit of the a squared unit of the plate used; or group
	# them by pii just for making it easier to the demand constraint.

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

# Include all the code implementing the pricings referred below.
include("./pricing.jl")

function _no_arg_check_build_model(
	model, d :: Vector{D}, p :: Vector{P}, l :: Vector{S}, w :: Vector{S},
	L :: S, W :: S, options :: Dict{String, Any} = Dict{String, Any}()
) :: ByproductPPG2KP{D, S, P} where {D, S, P}
	bp = build_complete_model(model, d, p, l, w, L, W, options)
	debug = options["verbose"] && !options["quiet"]
	pricing = options["pricing"]
	debug && begin
		length_pe_before_pricing = length(model[:picuts])
		length_cm_before_pricing = length(model[:cuts_made])
		@show length_pe_before_pricing
		@show length_cm_before_pricing
	end
	if pricing != "none"
		@timeit "pricing" begin
			if pricing == "expected"
				pricing = (options["faithful2furini2016"] ? "furini" : "becker")
			end
			if pricing == "furini"
				bp = _furini_pricing!(
					model, bp, d, p, options["pricing-heuristic-seed"],
					options["pricing-alpha"], options["pricing-beta"], debug
				)
			else
				@assert pricing == "becker"
				bp = _becker_pricing!(
					bp, model, d, p, options["pricing-heuristic-seed"], debug
				)
			end
		end
	end
	if debug
		println("length_pe_after_pricing = $(length(model[:picuts]))")
		println("length_cm_after_pricing = $(length(model[:cuts_made]))")
	end
	if !options["no-unreachable-check"]
		@timeit "unreachable" bp = _purge_unreachable!(bp, model, debug)
	end

	return bp
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

