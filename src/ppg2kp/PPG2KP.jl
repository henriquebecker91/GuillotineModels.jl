module PPG2KP
# Include submodules. Even if not used here, they need to be here to be
# avaliable for users to access/import.
include("./Heuristic.jl")
include("./Args.jl")
include("./Enumeration.jl")

using UnPack # for @unpack
using DocStringExtensions # For better struct documentation
import ..TIMER # Global module timer for use with TimerOutputs.@timeit.
import ..throw_if_timeout_now
using .Enumeration # Where all the plate enumeration logic is defined.
export ByproductPPG2KP # re-export ByproductPPG2KP from Enumeration

using ..Utilities
using RandomNumbers.Xorshifts # for Xoroshiro128Plus (for heuristic rng)
import ..to_pretty_str # for debug

# Import third-party libraries.
using JuMP
import TimerOutputs.@timeit
import MathOptInterface
const MOI = MathOptInterface

@enum Dimension begin
	LENGTH = 1
	WIDTH = 2
end
Base.getindex(tuple :: Tuple, s :: Dimension) = tuple[Int(s)]

@enum Relation begin
	PARENT = 1
	FIRST_CHILD = 2
	SECOND_CHILD = 3
end
Base.getindex(tuple :: Tuple, s :: Relation) = tuple[Int(s)]

@enum BaseModel BECKER FURINI

# Include code that is focused on implementing a specific feature, but does
# not merit a module of their own. Non-exported methods of these files are not
# expected to be used in this module, but exported methods are exports from
# this module (ou specialization of parent modules).
include("./get_cut_pattern.jl")
include("./MIP_start.jl")

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
	# Update the variable names to make the debug to make more sense.
	for (i, var) in enumerate(to_keep_vars)
		set_name(var, "$(var_set_name_in_model)[$(i)]")
	end
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

"""
Auxiliar struct used to build the PPG2KP model.

Aggregates the inverted indexes that indicate the variables with non-zero
coefficients in model constraints.

$TYPEDFIELDS

"""
struct VarInvIndexes{P}
# Note the field order, as all of them have the same type but different
# meanings, switching the order in the default constructor can be calamitous.
"""
If indexed by a piece index, returns the indexes of all extraction variables
(picuts/np) that extract that piece from some plate.
"""
pii2pair :: Vector{Vector{P}}
"""
If indexed by a plate index, returns the indexes of all extraction variables
(picuts/np) that extract some piece from that plate.
"""
pli2pair :: Vector{Vector{P}}
"""
If indexed by a plate index, returns the indexes of all cut variables
(cuts_made/cuts) that have that plate as either first or second child (i.e.,
one of the plates generated by the cut).
"""
child2cut :: Vector{Vector{P}}
"""
If indexed by a plate index, returns the indexes of all cut variables
(cuts_made/cuts) that have that plate as the parent of the cut (i.e., the plate
that is being cut).
"""
parent2cut :: Vector{Vector{P}}
end

"""
    VarInvIndexes(num_piece_types, num_plate_types, np, cuts)

Utility constructor that creates the struct from a ByproductPPG2KP.
"""
@timeit TIMER function VarInvIndexes(
	num_piece_types :: P, num_plate_types :: P,
	@nospecialize(np), @nospecialize(cuts)
) :: VarInvIndexes{P} where {P <: Integer}
	pii2pair = [Vector{P}() for _ = 1:num_piece_types]
	pli2pair = [Vector{P}() for _ = 1:num_plate_types]
	child2cut = [Vector{P}() for _ = 1:num_plate_types]
	parent2cut = [Vector{P}() for _ = 1:num_plate_types]

	# Initialize all inverse indexes.
	for i in eachindex(cuts)
		parent, fchild, schild = cuts[i]
		push!(parent2cut[parent], i)
		push!(child2cut[fchild], i)
		!iszero(schild) && push!(child2cut[schild], i)
	end

	for i in eachindex(np)
		pli, pii = np[i]
		push!(pli2pair[pli], i)
		push!(pii2pair[pii], i)
	end

	# Note the parameter order, as all of them have the same type but different
	# meanings, switching the order can be calamitous.
	return VarInvIndexes(pii2pair, pli2pair, child2cut, parent2cut)
end

"""
    VarInvIndexes(bp :: ByproductPPG2KP{D, S, P}) :: VarInvIndexes{P}

Utility constructor that creates the struct from a ByproductPPG2KP.
"""
function VarInvIndexes(
	bp :: ByproductPPG2KP{D, S, P}
) :: VarInvIndexes{P} where {D, S, P}
	num_piece_types = convert(P, length(bp.l))
	num_plate_types = convert(P, length(bp.pli_lwb))
	return VarInvIndexes(num_piece_types, num_plate_types, bp.np, bp.cuts)
end

# TODO: change enumeration to also use a Dict? and then use
# something like create_normalized_arg_subset in a subset of PPG2KP args
@timeit TIMER function _gen_cuts_wo(
	::Type{P}, d :: Vector{D}, l :: Vector{S}, w :: Vector{S},
	L :: S, W :: S, options :: Dict{String, Any} = Dict{String, Any}()
) where {D, S, P}
	return gen_cuts(P, d, SortedLinkedLW(D, l, w), L, W;
		ignore_2th_dim = options["ignore-2th-dim"],
		ignore_d = options["ignore-d"],
		round2disc = options["round2disc"],
		no_cut_position = options["no-cut-position"],
		no_redundant_cut = options["no-redundant-cut"],
		no_furini_symmbreak = options["no-furini-symmbreak"],
		faithful2furini2016 = options["faithful2furini2016"],
		quiet = options["quiet"], verbose = options["verbose"]
	)
end

@timeit TIMER function _build_base_model!(
	model, d :: Vector{D}, p :: Vector{P},
	bp :: ByproductPPG2KP{D, S, P}, inv_idxs :: VarInvIndexes{P},
	options :: Dict{String, Any} = Dict{String, Any}()
) where {D, S, P}
	# shorten the names
	@unpack l, w, L, W, np, cuts, pli_lwb = bp
	@unpack pii2pair, pli2pair, child2cut, parent2cut = inv_idxs
	num_piece_types = convert(D, length(d))
	num_plate_types = length(pli_lwb)

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
		0 <= cuts_made[i = 1:length(cuts)] <= pli_lwb[cuts[i][PARENT]][3]
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
    build_complete_model(model, d, p, l, w, L, W[, start]; options)

TODO: describe method.
Note: `model[:plate_cons]` is a container of constraints of type and set
`(GenericAffExpr{Float64,VariableRef}, LessThan{Float64})` that represent all
the plate types (i.e., it will be the same length as ByproductPPG2KP.pli_lwb
and each position in `model[:plate_cons]` will correspond to the plate
described in the same position of `ByproductPPG2KP.pli_lwb`).
"""
function build_complete_model(
	model, d :: Vector{D}, p :: Vector{P}, l :: Vector{S}, w :: Vector{S},
	L :: S, W :: S, start :: Float64 = time(),
	options :: Dict{String, Any} = Dict{String, Any}()
) :: ByproductPPG2KP{D, S, P} where {D, S, P}
	bp = _gen_cuts_wo(P, d, l, w, L, W, options)
	build_complete_model(model, d, p, bp, start, options)
	return bp
end

function build_complete_model(
	model, d :: Vector{D}, p :: Vector{P}, bp :: ByproductPPG2KP{D, S, P},
	start :: Float64 = time(), options :: Dict{String, Any} = Dict{String, Any}()
) :: ByproductPPG2KP{D, S, P} where {D, S, P}
	inv_idxs = VarInvIndexes(bp)
	limit :: Float64 = options["building-time-limit"]
	throw_if_timeout_now(start, limit)
	@assert issorted(
		@view bp.cuts[1:(bp.first_vertical_cut_idx-1)];
		by = c -> c[PARENT]
	)
	@assert	issorted(
		@view bp.cuts[bp.first_vertical_cut_idx:end];
		by = c -> c[PARENT]
	)
	_build_base_model!(model, d, p, bp, inv_idxs, options)

	return bp
end

# Include all the code implementing the pricings referred below.
include("./pricing.jl")

function _no_arg_check_build_model(
	model, d :: Vector{D}, p :: Vector{P}, l :: Vector{S}, w :: Vector{S},
	L :: S, W :: S, options :: Dict{String, Any}
) :: ByproductPPG2KP{D, S, P} where {D, S, P}
	# Define 'start', this is the point where the building-time-limit
	# starts to be counted.
	start :: Float64 = time()
	limit :: Float64 = options["building-time-limit"]
	# Enumerate plates and cuts.
	bp = _gen_cuts_wo(P, d, l, w, L, W, options)
	throw_if_timeout_now(start, limit)

	# Deal with pricing or its absence. Builds the model based on the enumerated
	# plates and cuts (or a subset of them).
	verbose = options["verbose"] && !options["quiet"]
	switch_method = options["Gurobi-LP-method-inside-furini-pricing"]
	pricing = options["pricing"]
	if pricing == "expected"
		pricing = (options["faithful2furini2016"] ? "furini" : "becker")
	end
	# The warnings are done here to avoid putting them inside a method
	# that will not be called anyway.
	current_solver_name = solver_name(model)
	if !options["quiet"] && switch_method > -2
		if pricing != "furini"
			@warn "The flag Gurobi-LP-method-inside-furini-pricing has a" *
				" non-default value, but the pricing is '$pricing'. This flag" *
				" will have no effect besides warnings."
		end
		if current_solver_name != "Gurobi"
			@warn "The flag Gurobi-LP-method-inside-furini-pricing has a" *
				" non-default value, but the solver is '$current_solver_name'." *
				" This flag will have no effect besides warnings."
		end
	end
	if verbose
		println("length_pe_before_pricing = $(length(bp.np))")
		println("length_cm_before_pricing = $(length(bp.cuts))")
		println("length_pc_before_pricing = $(length(bp.pli_lwb))")
	end
	bm = (options["faithful2furini2016"] ? FURINI : BECKER) :: BaseModel
	if pricing != "none"
		# It is important to note that the pricing and the bm (BaseModel) have the
		# same two values ("becker"/BECKER and "furini"/FURINI) but they are
		# orthogonal/independent. You can have a FURINI model with "becker" pricing
		# and vice-versa. It is just because each author has defined both a model
		# and an optional pricing method (to be executed over it) but both pricings
		# are compatible with both models, and we allow them to be mixed. The
		# pricing methods just need to know the underlying base model to be able to
		# MIP-start it corretly (initially at least, this comment may be out of
		# date).
		if pricing == "furini"
			bp = _furini_pricing!(model, bp, d, p, start, options)
			verbose && (pricing_time = past_section_seconds(TIMER, "_furini_pricing!"))
		else
			@assert pricing == "becker"
			bp = _becker_pricing!(model, bp, d, p, start, options)
			verbose && (pricing_time = past_section_seconds(TIMER, "_becker_pricing!"))
		end
		verbose && @show pricing_time
	else # in the case there is no pricing phase
		# Needs to build the mode here, as when there is pricing the pricing
		# procedure is responsible for building the model.
		build_complete_model(model, d, p, bp, start, options)
		if mip_start == "guaranteed" # we have to do the MIP-start ourselves
			heuristic_seed = options["heuristic-seed"]
			mip_start_by_heuristic!(model, bp, d, p, heuristic_seed, bm)
		end
	end
	if verbose
		println("length_pe_after_pricing = $(length(model[:picuts]))")
		println("length_cm_after_pricing = $(length(model[:cuts_made]))")
		println("length_pc_after_pricing = $(length(model[:plate_cons]))")
	end
	throw_if_timeout_now(start, limit) # after the pricings or warm-starts

	# Check if variables and constraints made useless by the deletion
	# of other variables will be kept or not.
	purge = !options["do-not-purge-unreachable"]
	bp = _handle_unreachable!(bp, model, verbose, purge)
	if verbose && purge
		println("length_pe_after_purge = $(length(model[:picuts]))")
		println("length_cm_after_purge = $(length(model[:cuts_made]))")
		println("length_pc_after_purge = $(length(model[:plate_cons]))")
	end
	throw_if_timeout_now(start, limit) # after handle_unreachable!

	return bp
end

# Only used to define `build_model(Val{:PPG2KP}, ...)`.
import ..Utilities.Args.create_normalized_arg_subset
import ..Utilities.Args.accepted_arg_list
import ..build_model
"""
    build_model(::Val{:PPG2KP}, model, d, p, l, w, L, W[, options])


"""
@timeit TIMER function build_model(
	::Val{:PPG2KP}, model, d :: Vector{D}, p :: Vector{P},
	l :: Vector{S}, w :: Vector{S}, L :: S, W :: S,
	options :: Dict{String, Any} = Dict{String, Any}()
) :: ByproductPPG2KP{D, S, P} where {D, S, P}
	norm_options = create_normalized_arg_subset(
		options, accepted_arg_list(Val(:PPG2KP))
	)
	bmr = _no_arg_check_build_model(model, d, p, l, w, L, W, norm_options)
	# show(stdout, "text/plain", bmr.cuts)
	# show(stdout, "text/plain", bmr.pli_lwb)
	return bmr
end

end # module

