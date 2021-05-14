module PPG2KP

# Constant used for the problems that need no change in the enumeration
# process (and often share the same method bodies).
const SIMILAR_4 = Union{Val{:G2KP}, Val{:G2MKP}, Val{:G2OPP}, Val{:G2CSP}}

# Include submodules. Even if not used here, they need to be here to be
# avaliable for users to access/import.
include("./Heuristic.jl")
include("./Args.jl")
include("./Enumeration.jl")

# Everything that is exported by PPG2KP.Enumeration ends up being used here.
using .Enumeration # Where all the plate enumeration logic is defined.
export ByproductPPG2KP # re-export ByproductPPG2KP from Enumeration

# Import third-party libraries.
using JuMP
using UnPack # for @unpack
import TimerOutputs.@timeit
import MathOptInterface
using RandomNumbers.Xorshifts # for Xoroshiro128Plus (for heuristic rng)
const MOI = MathOptInterface
using DocStringExtensions # For better struct documentation

# Import globals/modules/functions defined in GuillotineModels (main module).
import ..TIMER # Global module timer for use with TimerOutputs.@timeit.
import ..throw_if_timeout_now
import ..to_pretty_str # for debug
import ..CutPattern # type returned by get_cut_pattern
import ..BuildStopReason, ..BUILT_MODEL, ..FOUND_OPTIMUM
using ..Utilities
import ..Data.G2KP
import ..Data.SLOPP
import ..Data.MHLOPPW
import ..Data.SSSCSP
using ..Utilities: expand # for get_cut_pattern

# Included at the start because ModelByproduct has a RotationAwareData field.
include("rotation.jl")

struct ModelByproduct{D, S, P}
	preprocess_byproduct :: ByproductPPG2KP{D, S, P}
	found_optimum :: Bool
	# TODO: `optimum_if_found` will need to become a `Vector` of
	# `CutPattern`s in the future if this feature is to be supported by
	# G2CSP and G2MKP.
	optimum_if_found :: CutPattern{D, S}
	rad :: Union{Nothing, RotationAwareData{D, S}}
end

# Convenience in case of FOUND_OPTIMUM. `found_optimum` is set
# to `true` and `optimum_if_found` stores the solution.
function ModelByproduct(
	bp :: ByproductPPG2KP{D, S, P}, optimum :: CutPattern{D, S}
) where {D, S, P}
	return ModelByproduct{D, S, P}(bp, true, optimum)
end

# Convenience in case of BUILT_MODEL. `found_optimum` is set
# to `false` and `optimum_if_found` stores an empty (but valid) solution.
function ModelByproduct(bp :: ByproductPPG2KP{D, S, P}) where {D, S, P}
	return ModelByproduct{D, S, P}(bp, false, CutPattern(bp.L, bp.W), nothing)
end
# Convenience in case of BUILT_MODEL and rotation was allowed.
function ModelByproduct(
	bp :: ByproductPPG2KP{D, S, P}, rad :: RotationAwareData{D, S}
) where {D, S, P}
	return ModelByproduct{D, S, P}(bp, false, CutPattern(bp.L, bp.W), rad)
end

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
		bp.cuts, new_fvci, bp.np, bp.pli_lwb, bp.d, bp.l, bp.w, bp.L, bp.W
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
		!iszero(fchild) && push!(child2cut[fchild], i)
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

#=
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
=#

@timeit TIMER function _gen_cuts_wo(
	::Type{P}, d :: Vector{D}, sllw :: SortedLinkedLW{D, S},
	L :: S, W :: S, options :: Dict{String, Any} = Dict{String, Any}()
) where {D, S, P}
	return gen_cuts(P, d, sllw, L, W,
		Val(options["allow-rotation"] & options["mirror-plates"]);
		ignore_2th_dim = options["ignore-2th-dim"],
		ignore_d = options["ignore-d"],
		round2disc = options["round2disc"],
		hybridize_with_restricted = options["hybridize-with-restricted"],
		aggressive_hybridization = options["aggressive-hybridization"],
		no_cut_position = options["no-cut-position"],
		no_redundant_cut = options["no-redundant-cut"],
		no_furini_symmbreak = options["no-furini-symmbreak"],
		faithful2furini2016 = options["faithful2furini2016"],
		quiet = options["quiet"], verbose = options["verbose"],
	)
end

# INTERNAL USE.
#
# (Only used in _build_base_model! to make the shared demand constraint.)
#
# Helper functions that make it easier to index collections using the
# Vector{Union{D, Tuple{D, D}}} index lists.
function _flatten_index(a, i :: Integer, f)
	return a[i]
end
function _flatten_index(a, i :: Tuple{I, I}, f) where {I <: Integer}
	return f(a[first(i)], a[last(i)])
end

# INTERNAL USE.
#
# (Only used in _build_base_model! to compute a cut upper bound.)
#
# Helper function that computes a loose upper bound on the number of original
# plates available.
_ub_plates(inst :: MHLOPPW) = only(inst.available)
function _ub_plates(inst :: SSSCSP{D, S, P}) where {D, S, P}
	@unpack L, W, l, w, d = inst
	# The simplest (but terribly loose) upper bound is just `sum(d)`
	# (i.e., there is a solution in which each piece is in a plate of its
	# own). The following is the second-simplest thing: we have 4 buckets,
	# for pieces in which li <= L÷8 and wi <= W÷8 (we can fit 64 of these
	# in a plate), in which li <= L÷4 and wi <= W÷4 (we can fit 16 of these
	# in a plate), in which li <= L÷2 and wi <= W÷2 (we can fit 4 of these
	# in a plate), and in which li > L÷2 or wi > W÷2 (we give a whole plate
	# for these). This way, the bound is not terrible in the presence of
	# small- and medium-sized pieces.
	le8qt64, le4qt16, le2qt4, gt2 = zero(D), zero(D), zero(D), zero(D)

	L8, L4, L2 = L ÷ 8, L ÷ 4, L ÷ 2
	W8, W4, W2 = W ÷ 8, W ÷ 4, W ÷ 2
	for (li, wi, di) in zip(l, w, d)
		if     li <= L8 && wi <= W8
			le8qt64 += di
		elseif li <= L4 && wi <= W4
			le4qt16 += di
		elseif li <= L2 && wi <= W2
			le2qt4  += di
		else # i.e., li > L2 || wi > W2
			gt2     += di
		end
	end

	# The trick below is: if we sum (divisor - 1) to the number about to be
	# divided we transform ÷ into a round-up (not round-down) integer division.
	final_ub = (
		((le8qt64 + 63) ÷ 64) +
		((le4qt16 + 15) ÷ 16) +
		((le2qt4  +  3) ÷  4) +
		gt2
	)

	@assert final_ub <= sum(d)

	return final_ub
end

@timeit TIMER function _build_base_model!(
	problem :: SIMILAR_4, inst, model, bp :: ByproductPPG2KP{D, S, P},
	inv_idxs :: VarInvIndexes{P},
	rad :: Union{Nothing, RotationAwareData{D, S}},
	options :: Dict{String, Any} = Dict{String, Any}();
	build_LP_not_MIP = false
	# `build_LP_not_MIP` is used internally by the pricing (the relaxation needs
	# to be solved to solve the problem as a MIP).
	# `options["relax2lp"]` is a generic option to be handled in the generic
	# (i.e., formulation-agnostic) part of the code.
) where {D, S, P}
	# shorten the names
	@unpack d, l, w, L, W, np, cuts, cut_extraction, pli_lwb = bp
	@unpack pii2pair, pli2pair, child2cut, parent2cut = inv_idxs
	num_piece_types = convert(D, length(d))
	num_plate_types = length(pli_lwb)

	# This probably should not be here (but in an outer function)
	# but for now this will suffice.
	pii2dci = Vector{D}[Vector{D}() for i in 1:num_piece_types]
	for (dci, pii) in pairs(cut_extraction)
		!iszero(pii) && push!(pii2dci[pii], dci)
	end
	#@show length(cuts)
	#@show length(cut_extraction)

	if problem == Val(:G2KP) || problem == Val(:G2MKP)
		if options["allow-rotation"]
				# If a rotation-aware piece (i.e., dummy) maps to two original
				# pieces (only happens if both are perfect rotations of each other)
				# then we can map the dummy index to the any of the original pieces
				# (here we choose the first) because both share the same profit value
				# (this is a requirement to be a perfect rotation).
				rai2opi = first.(rad.sdi2opi[rad.rai2sdi]) :: Vector{D}
				p = inst.p[rai2opi] :: Vector{P}
		else
				p = inst.p :: Vector{P}
		end
	end

	# Note: using `d` as upper bound of `picuts` and `dc_sells` works even
	# with rotation enabled because `d` is the `rad.shared_d` expanded for the
	# dummies.
	@variable(model,
		0 <= picuts[i = 1:length(np)] <= d[np[i][2]], integer = !build_LP_not_MIP
	)
	if options["hybridize-with-restricted"]
		# dc_sells: Means "double-cut sells". For each piece type `pii` (or dummy
		# piece type in the case of rotation enabled), we have `dc_sells[pii]`
		# storing the number of pieces obtained by means of double cuts that were
		# actually sold. This variable is necessary because we cannot assume the
		# piece obtained by a double cut is always sold.  When double cuts are
		# enabled (i.e., hybridize-with-restricted), the normal cuts replaced by
		# double cuts cease to exist, and the double cuts may be needed in some
		# situations just to trim waste, without selling the piece they extract
		# because that piece demand is already exhausted.
		@variable(model,
			0 <= dc_sells[i = 1:num_piece_types] <= d[i],
			integer = !build_LP_not_MIP
		)
	end

	# We prefer to leave `cuts_made` without an upper bound to risk
	# making a bound that will not work for some combination of parameters.
	# EDIT: we are adding back a loose upper bound because it seems like
	# Gurobi have a bug that may give wrong answers if upper bounds are not set.
	# The `ub_cut_in_single_plate` is just an upper bound on the amount of
	# plates p that fit in a single original plate, where p is the parent plate
	# of the cut (i.e., the plate that is being divided by the cut).
	ub_cut_in_single_plate = getindex.(pli_lwb, 3)[getindex.(cuts, 1)]
	ub_cut = if problem === Val(:G2MKP) || problem === Val(:G2CSP)
		# NOTE: _ub_plates assumes no-rotation because this is a valid upper bound
		# for both rotation and no-rotation variants.
		ub_cut_in_single_plate .* _ub_plates(inst)
	else
		ub_cut_in_single_plate
	end
	@variable(model,
		0 <= cuts_made[i = 1:length(cuts)] <= ub_cut[i],
		integer = !build_LP_not_MIP
	)

	if options["hybridize-with-restricted"]
		# Constraint that guarantee dc_sells only sell pieces actually extracted
		# by double cuts.
		@constraint(model,
			dc_link_con[i = 1:length(d)], dc_sells[i] <= sum(cuts_made[pii2dci[i]])
		)
	end

	# If the problem is a knapsack problem, then the objective is to maximize
	# the profit from the packed/extracted pieces.
	# If the problem is CSP, then the objective is to minimize the number of
	# original plates used (all pieces must be packed).
	# If it is OPP, then there is no objective function.
	if problem === Val(:G2KP) || problem === Val(:G2MKP)
		# If the objective is profit maximization we need to account for the
		# sold pieces obtained by double cut extractionsa (`dc_sells`).
		if options["hybridize-with-restricted"]
			@objective(model, Max,
				sum(p[pii] * sum(picuts[pii2pair[pii]]) for pii = 1:num_piece_types) +
				sum(p[pii] * dc_sells[pii] for pii = 1:num_piece_types)
			)
		else
			@objective(model, Max,
				sum(p[pii] * sum(picuts[pii2pair[pii]]) for pii = 1:num_piece_types)
			)
		end
	elseif problem === Val(:G2CSP)
		# We create a variable to represent the number of original plates used:
		# `b`. This variable is not strictly necessary, we could instead minimize
		# the cuts and extractions made over plates of plate type #1. However,
		# without the variable we will have one less constraint compared to the
		# other problems (there will not be the `plate_cons` #1 constraint). This
		# would add some unfortunate assymetry; either the list of constraint would
		# start at the second index, or to map a plate to its constraint we would
		# need to always subtract one of the plate index (when it is G2CSP,
		# otherwise not, therefore a mess). This seems to be the better design,
		# no code (even G2CSP-specific) really need to be aware of this variable.
		@variable(model, b, integer = !build_LP_not_MIP)
		#@objective(model, Min, sum([b])) # the sum prevents errors with GLPK
		@objective(model, Min, b)
	end

	# c1: # Either there is a limited amount of original plates available
	# (one if it is KP or OPP, and 'n' if it is MKP), or there is an
	# unlimited amount of original plates available (CSP) but we are minimizing
	# the amount actually used (variable `b`).
	#
	# The original plates may be used either to extract pieces directly
	# or to be cut into two new subplates.
	if problem === Val(:G2CSP)
		plate_number_one = @constraint(model,
			sum(picuts[pli2pair[1]]) + sum(cuts_made[parent2cut[1]]) <= b
		)
	else
		qt_original :: Int = problem === Val(:G2MKP) ? only(inst.available) : 1
		plate_number_one = @constraint(model,
			sum(picuts[pli2pair[1]]) + sum(cuts_made[parent2cut[1]]) <= qt_original
		)
	end

	# c2: for each subplate type that is not the original plate, such subplate
	# type will be available the number of times it was the child of a cut,
	# subtracted the number of times it had a piece extracted or used for
	# further cutting. This constraint is the basis of the model and similar
	# for all problems it supports.
	# NOTE: the range is 1:(num_plate_types-1) and not 2:num_plate_types
	# because if the range does not start at 1, the returned type is not
	# a vector but a special JuMP type.
	plate_numbers_2_to_m = @constraint(model, [ppli=1:(num_plate_types-1)],
		sum(picuts[pli2pair[ppli+1]]) + sum(cuts_made[parent2cut[ppli+1]]) <=
		sum(cuts_made[child2cut[ppli+1]])
	)
	# give the constraints a name
	set_name.(plate_numbers_2_to_m, "plate_cons" .* string.(2:num_plate_types))
	# put them in the expected key inside the model
	model[:plate_cons] = pushfirst!(plate_numbers_2_to_m, plate_number_one)

	# The demand constraint, yet simple, is modified in complementary ways by
	# basically every flag that changes the formulation, leading an exponential
	# explosion of slightly different demand constraints. To avoid 3+ levels of
	# nested `if`s, we try to build each part of the demand constraint
	# separatedly (depending only in the options that affect that single part)
	# and then combine the parts at the end.

	# The right-hand-side of the demand constraint is determined by the
	# (no-)rotation. The (no-)rotation also defines how we translate
	# demand constraint indexes to piece indexes (`di2pii`).
	if options["allow-rotation"]
		demand_con_rhs = rad.shared_demand
		di2pii = rad.sdi2rai :: fieldtype(RotationAwareData{D, S}, :sdi2rai)
	else
		demand_con_rhs = d
		# We use typeof to avoid type-unstability on binding `di2pii`.
		di2pii = fieldtype(RotationAwareData{D, S}, :sdi2rai)()
		append!(di2pii, eachindex(d))
	end
	# With the rhs and the translator we may build lhs expression outside
	# of any `if`.
	demand_con_lhs = @expression(model, [di=1:length(di2pii)],
		sum(picuts[_flatten_index(pii2pair, di2pii[di], vcat)])
	)
	# If we are hybridizing with restricted we need to update lhs to
	# also include the extractions/sells made by double cuts.
	if options["hybridize-with-restricted"]
		demand_con_lhs = @expression(model, [di=1:length(di2pii)],
			demand_con_lhs[di] + sum(dc_sells[collect(di2pii[di])])
		)
	end
	# Now, with lhs and rhs in hand, we will build contraints that only differ in
	# the sense (i.e., either `<=` or `>=`).
	# For knapsack problems, demand means an upper bound on how many of each
	# piece type may be extracted/packed, for OPP and CSP it means a lower bound
	# on the number of each piece type that *must* be extracted/packed.
	if problem === Val(:G2KP) || problem === Val(:G2MKP)
		@constraint(model,
			demand_con[di=1:length(di2pii)],
			demand_con_lhs[di] <= demand_con_rhs[di]
		)
	elseif problem === Val(:G2OPP) || problem === Val(:G2CSP)
		@constraint(model,
			demand_con[di=1:length(di2pii)],
			demand_con_lhs[di] >= demand_con_rhs[di]
		)
	end

	#= Disabled while we rework it in a problem-agnostic fashion.
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
	=#

	return
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
#=
function build_complete_model(
	model, d :: Vector{D}, p :: Vector{P}, l :: Vector{S}, w :: Vector{S},
	L :: S, W :: S, start :: Float64 = time(),
	options :: Dict{String, Any} = Dict{String, Any}()
) :: ByproductPPG2KP{D, S, P} where {D, S, P}
	bp = _gen_cuts_wo(P, d, l, w, L, W, options)
	build_complete_model(model, p, bp, start, options)
	return bp
end
=#
function build_complete_model(
	problem, instance, model, bp :: ByproductPPG2KP{D, S, P},
	rad :: Union{Nothing, RotationAwareData{D, S}},
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
	_build_base_model!(problem, instance, model, bp, inv_idxs, rad, options)

	return bp
end

# Include all the code implementing the pricings referred below.
include("./pricing.jl")

# Auxiliar function to get the geometric info from the homogeneous
# problems and pass it to the enumeration.
function _homo_geometry(
	i :: Union{G2KP{D, S, P}, SSSCSP{D, S, P}}
) where {D, S, P}
	return D, S, P, i.d, i.l, i.w, i.L, i.W
end
function _homo_geometry(i :: MHLOPPW{D, S, P}) where {D, S, P}
	return D, S, P, i.dub, i.l, i.w, only(i.L), only(i.W)
end
function _get_profit(
	i :: Union{G2KP{D, S, P}, MHLOPPW{D, S, P}}
) :: Vector{P} where {D, S, P}
	return i.p
end
function _get_profit(
	i :: SSSCSP{D, S, P}
) :: Vector{P} where {D, S, P}
	return zeros(P, length(i.d))
end

function _no_arg_check_build_model(
	problem, instance, model, options :: Dict{String, Any}
) :: Tuple{BuildStopReason, <: Any}
	# Define 'start', this is the point where the building-time-limit
	# starts to be counted.
	start :: Float64 = time()
	limit :: Float64 = options["building-time-limit"]
	verbose = options["verbose"] && !options["quiet"]
	# Enumerate plates and cuts.
	enumeration_time = @elapsed begin
		D, _, P, d, l, w, L, W = _homo_geometry(instance)
		sllw = SortedLinkedLW(D, l, w)
		if options["allow-rotation"]
			rad = build_RAD(L, W, sllw, d, _get_profit(instance))
			# Here we expand the shared demand to all dummy pieces it applies,
			# so it is invisible to _gen_cuts_wo that we are dealing with
			# a rotation variant.
			unshared_demand = rad.shared_demand[rad.rai2sdi]
			bp = _gen_cuts_wo(P, unshared_demand, rad.ra_sllw, L, W, options)
		else
			bp = _gen_cuts_wo(P, d, sllw, L, W, options)
		end
	end
	if verbose
		@show enumeration_time
		println("qt_pevars_after_preprocess = $(length(bp.np))")
		println("qt_cmvars_after_preprocess = $(length(bp.cuts))")
		println("qt_plates_after_preprocess = $(length(bp.pli_lwb))")
	end
	throw_if_timeout_now(start, limit)

	# Deal with pricing or its absence. Builds the model based on the enumerated
	# plates and cuts (or a subset of them).
	switch_method = options["Gurobi-LP-method-inside-furini-pricing"]
	pricing = options["pricing"]
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
	bm = (options["faithful2furini2016"] ? FURINI : BECKER) :: BaseModel
	if pricing != "none"
		problem !== Val(:G2KP) && @error(
			"The pricing was only implemented for the G2KP problem. Other" *
			" problems have no implementation of the pricing."
		)
		options["allow-rotation"] && @error(
			"The pricing code was not yet examined to check if it supports" *
			" rotation. The combination of any pricing method plus rotation" *
			" is disabled for now."
		)
		# It is important to note that the pricing and the bm (BaseModel) have the
		# same two values ("becker"/BECKER and "furini"/FURINI) but they are
		# orthogonal/independent. You can have a FURINI model with "becker" pricing
		# and vice-versa. It is just because each author has defined both a model
		# and an optional pricing method (to be executed over it) but both pricings
		# are compatible with both models, and we allow them to be mixed. The
		# pricing methods just need to know the underlying base model to be able to
		# MIP-start it corretly (initially at least, this comment may be out of
		# date).
		pricing_time = @elapsed begin
			if pricing == "furini"
				bsr, mbp = _furini_pricing!(model, bp, p, start, options)
				bp = mbp.preprocess_byproduct
			else
				@assert pricing == "becker"
				bp = _becker_pricing!(model, bp, p, start, options)
			end
		end
		verbose && @show pricing_time
		# Could be cleaner...
		if @isdefined(bsr) && bsr == FOUND_OPTIMUM
			@assert @isdefined(mbp)
			return bsr, mbp
		end
	else # in the case there is no pricing phase
		# Needs to build the mode here, as when there is pricing the pricing
		# procedure is responsible for building the model.
		build_complete_model(
			problem, instance, model, bp,
			(options["allow-rotation"] ? rad : nothing), start, options
		)
		if options["MIP-start"] == "guaranteed"
			problem !== Val(:G2KP) && @error(
				"The MIP-start was only implemented for the G2KP problem. Other" *
				" problems have no implementation of the MIP-start."
			)
			options["allow-rotation"] && @error(
				"The MIP-start was not yet examined to check if it supports" *
				" rotation. Using MIP-start plus rotation is disabled for now."
			)
			heuristic_lb_time = @elapsed begin
				(heuristic_lb, _, _), _ = mip_start_by_heuristic!(
					model, bp, p, options["heuristic-seed"], bm
				)
			end
			if verbose
				@show heuristic_lb_time
				@show heuristic_lb
			end
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
	# TODO: if we ever start dealing with multiple original plates, we may
	# need to update _handle_unreachable! (the exception is if we go with
	# some maximum.((L, W)) workaround were the original plates are all
	# cut from a dummy root plate).
	purge_unreachable_time = @elapsed begin
		purge = !options["do-not-purge-unreachable"]
		bp = _handle_unreachable!(bp, model, verbose, purge)
		if verbose && purge
			println("qt_pevars_after_purge = $(length(model[:picuts]))")
			println("qt_cmvars_after_purge = $(length(model[:cuts_made]))")
			println("qt_plates_after_purge = $(length(model[:plate_cons]))")
		end
	end
	verbose && @show purge_unreachable_time
	throw_if_timeout_now(start, limit) # after handle_unreachable!

	options["allow-rotation"] && return BUILT_MODEL, ModelByproduct(bp, rad)

	return BUILT_MODEL, ModelByproduct(bp)
end

# Only used to define `build_model(Val{:PPG2KP}, ...)`.
import ..Utilities.Args.create_normalized_arg_subset
import ..Utilities.Args.accepted_arg_list
import ..build_model
#=
"""
    build_model(::Val{:G2KP}, ::Val{:PPG2KP}, instance::G2KP, model[, options])
    build_model(::Val{:G2MKP}, ::Val{:PPG2KP}, instance::MHLOPPW, model[, options])
    build_model(::Val{:G2CSP}, ::Val{:PPG2KP}, instance::SSSCSP, model[, options])
    build_model(::Val{:G2OPP}, ::Val{:PPG2KP}, instance::SSSCSP, model[, options])

Build a PPG2KP-style model for a G2KP `instance` inside `model`.

Changes `model` by adding variables and constraints. `options` takes
arguments described in `Utilities.Args.accepted_arg_list(::Val{:PPG2KP})`.
"""
=#
#=
@timeit TIMER function build_model(
	HOMO_PROBLEMS, ::Val{:PPG2KP}, instance :: G2KP{D, S, P}, model,
	options :: Dict{String, Any} = Dict{String, Any}()
) :: Tuple{BuildStopReason, <: Any} where {D, S, P}
	norm_options = create_normalized_arg_subset(
		options, accepted_arg_list(Val(:PPG2KP))
	)
	return _no_arg_check_build_model(model, norm_options)
end
=#

const HOMO_PAIRS = [
	(:G2KP, :G2KP), (:G2CSP, :SSSCSP), (:G2OPP, :SSSCSP)
]

for (problem, instance) in HOMO_PAIRS
	eval(quote
		@timeit TIMER function build_model(
			::Val{$(QuoteNode(problem))}, ::Val{:PPG2KP}, instance :: $instance{D, S, P}, model,
			options :: Dict{String, Any} = Dict{String, Any}()
		) :: Tuple{BuildStopReason, <: Any} where {D, S, P}
			norm_options = create_normalized_arg_subset(
				options, accepted_arg_list(Val(:PPG2KP))
			)
			return _no_arg_check_build_model(
				Val($(QuoteNode(problem))), instance, model, norm_options
			)
		end
	end)
end

"""
    build_model(::Val{:G2KP}, ::Val{:PPG2KP}, instance::SLOPP, model[, options])

Build a PPG2KP-style model for a G2KP `instance` inside `model`.

Changes `model` by adding variables and constraints. `options` takes
arguments described in `Utilities.Args.accepted_arg_list(::Val{:PPG2KP})`.

This is a convenience method that takes a SLOPP instance. It errors if
the `SLOPP` object has a non-zero value in the `dlb` field. It is
the same as calling with a `G2KP` instance but with the `dub` field
being used as the `d` field (of `G2KP`).
"""
@timeit TIMER function build_model(
	::Val{:G2KP}, ::Val{:PPG2KP}, instance :: SLOPP{D, S, P}, m,
	options :: Dict{String, Any} = Dict{String, Any}()
) :: Tuple{BuildStopReason, <: Any} where {D, S, P}
	@assert all(iszero, instance.dlb)
	@unpack p, l, w, L, W = instance
	d = instance.dub
	return build_model(
		Val(:G2KP), Val(:PPG2KP), G2KP(L, W, l, w, d, p), m, options
	)
end

"""
    build_model(::Val{:G2CSP}, ::Val{:PPG2KP}, instance::G2KP, model[, options])

Build a PPG2KP-style model for a G2CSP `instance` inside `model`.

Changes `model` by adding variables and constraints. `options` takes
arguments described in `Utilities.Args.accepted_arg_list(::Val{:PPG2KP})`.

This is a convenience method that takes a G2KP instance. It allows
easily testing G2CSP over the same instances used commonly for
Knapsack Problem (G2KP). If `options["quiet"]` is not passed,
this method will always show a warnig, as the instance type is
not ideal.
"""
@timeit TIMER function build_model(
	::Val{:G2CSP}, ::Val{:PPG2KP}, instance :: G2KP{D, S, P}, m,
	options :: Dict{String, Any} = Dict{String, Any}()
) :: Tuple{BuildStopReason, <: Any} where {D, S, P}
	@unpack L, W, l, w, d = instance
	if !options["quiet"]
		@warn "Attention: converting a G2KP instance to a SSSCSP instance to" *
			" solve it for the G2CSP (problem). Profits dropped."
	end
	return build_model(
		Val(:G2CSP), Val(:PPG2KP), SSSCSP{D, S, P}(L, W, l, w, d), m, options
	)
end

"""
    build_model(::Val{:G2KP}, ::Val{:PPG2KP}, instance::SSSCSP, model[, options])

Build a PPG2KP-style model for a G2OPP `instance` inside `model`.

Changes `model` by adding variables and constraints. `options` takes
arguments described in `Utilities.Args.accepted_arg_list(::Val{:PPG2KP})`.

THIS IS A HACK THAT SHOULD BE DELETED AFTER. It gets a SSSCSP instance, and
does not interpret it as a G2CSP/G2BPP instance, but instead as a G2OPP
instance and then it transforms it in a G2KP problem with unitaty profits,
so the G2OPP (a decision problem) is true if obj matches the sum(d), and
is false if obj is smaller than sum(d) (a value larger than sum(d) is
impossible).
"""
@timeit TIMER function build_model(
	::Val{:G2KP}, ::Val{:PPG2KP}, instance :: SSSCSP{D, S, P}, m,
	options :: Dict{String, Any} = Dict{String, Any}()
) :: Tuple{BuildStopReason, <: Any} where {D, S, P}
	@unpack L, W, l, w, d = instance
	if !options["quiet"]
		@warn "Attention: converting a SSSCSP (understood as a G2OPP) instance to a G2KP instance to" *
			" solve it for the G2KP (problem). Unitary profits added."
	end
	return build_model(
		Val(:G2KP), Val(:PPG2KP), G2KP{D, S, P}(L, W, l, w, d, ones(P, length(d))), m, options
	)
end

"""
    build_model(::Val{:G2MKP}, ::Val{:PPG2KP}, instance::MHLOPPW, model[, options])

Build a PPG2KP-style model for a G2MKP `instance` inside `model`.

Changes `model` by adding variables and constraints. `options` takes
arguments described in `Utilities.Args.accepted_arg_list(::Val{:PPG2KP})`.

For now, this is the only format accepted by G2MKP, even if the G2MKP
refers to the homogeneous variant and, therefore, does not accept
large objects of different dimensions (just multiple copies of the
same large objects). A method accepting a more adequate format
(i.e., that only allow a single type of large object) could be
implemented.
"""
@timeit TIMER function build_model(
	::Val{:G2MKP}, ::Val{:PPG2KP}, instance :: MHLOPPW{D, S, P}, m,
	options :: Dict{String, Any} = Dict{String, Any}()
) :: Tuple{BuildStopReason, <: Any} where {D, S, P}
	if (length(instance.L), length(instance.W)) != (1, 1)
		error(
			"The G2MKP problem refers to the homogeneous variant." *
			" Multiple original plates (i.e., large objects) of distinct" *
			" dimensions are not allowed."
		)
	end
	norm_options = create_normalized_arg_subset(
		options, accepted_arg_list(Val(:PPG2KP))
	)
	return _no_arg_check_build_model(
		Val(:G2MKP), instance, m, norm_options
	)
end

end # module

