module Utilities

using JuMP
using MathOptInterface
const MOI = MathOptInterface
using MathOptFormat
using TimerOutputs

# JuMP/'Mathematical Model' related utilities
export num_all_constraints, reduced_cost, delete_vars_by_pricing!
export relax_vars!, relax_all_vars!, unrelax_vars!, unrelax_all_vars!
export SavedBound, save_bound_if_exists!, restore_bound!
export which_vars_to_keep, fix_vars!
export save_model

# Style guideline: as the module block is left unindented, the @timeit
# blocks that wrap the whole method body also are not indented.

"""
    save_model(model, filename = "save_model.mps") :: Nothing

Save `model` using the MPS format in `filename`.
"""
function save_model(model, filename = "saved_model.mps") :: Nothing
	@timeit "save_model" begin
	mps_model = MathOptFormat.MPS.Model()
	MOI.copy_to(mps_model, backend(model))
	MOI.write_to_file(mps_model)
	end # timeit
	nothing
end

"""
    num_all_constraints(m) :: Int64

JuMP only allow to query the number of constraints of some specific type;
this method queries all constraint types used in the model and then sums the
number of constraints of each type.
"""
function num_all_constraints(m) :: Int64
	sum = 0 :: Int64
	for (ftype, stype) in list_of_constraint_types(m)
		sum += num_constraints(m, ftype, stype)
	end
	return sum
end

# see https://github.com/JuliaOpt/MathOptInterface.jl/issues/776
function reduced_cost(var) :: Float64
	rc = 0.0
	has_upper_bound(var) && (rc += shadow_price(UpperBoundRef(var)))
	has_lower_bound(var) && (rc += shadow_price(LowerBoundRef(var)))
	is_fixed(var) && (rc += shadow_price(FixRef(var)))
	!has_upper_bound(var) && !has_lower_bound(var) && !is_fixed(var) &&
		@warn "CAUTION: reduce_cost was called over a variable with no bounds"
	rc
end

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
# be extracted).
# TODO: check if this is necessary, and if it is, it should not go there,
# but inside PPG2KP module instead.
function which_vars_to_keep(model, lb :: Float64)
	obj = objective_value(model)
	rcs = reduced_cost.(all_variables(model))
	# We keep the ones equal to lb (not strictly necessary) to guarantee the
	# warm start will work (the variables needed for it should be there).
	return map(rc -> obj - rc >= lb, rcs)
end

"Stores a variable and its old bounds so they may be restored."
struct SavedBound
	"A reference to the variable."
	var :: VariableRef
	"Stores wether the variable was fixed."
	was_fixed :: Bool
	"If the variable was fixed, to which value they were fixed."
	fix_value :: Float64
	"Stores wether the variable had a lower bound."
	had_lb :: Bool
	"If the variable had a lower bound, the value of their lower bound."
	lb :: Float64
	"Stores wether the variable had an upper bound."
	had_ub :: Bool
	"If the variable had an upper bound, the value of their upper bound."
	ub :: Float64
end

"""
    restore_bound!(b :: SavedBound) :: Nothing

If `b.var` bounds are different than the ones specified by the rest
of the `b` fields, then change `b` to adhere to them.
"""
function restore_bound!(b :: SavedBound) :: Nothing
	if b.was_fixed && (!is_fixed(b.var) || fix_value(b.var) != b.fix_value)
		fix(b.var, b.fix_value; force = true)
	else
		is_fixed(b.var) && unfix(b.var)
		if b.had_lb && (!has_lower_bound(b.var) || lower_bound(b.var) != b.lb)
			set_lower_bound(b.var, b.lb)
		end
		if b.had_ub && (!has_upper_bound(b.var) || upper_bound(b.var) != b.ub)
			set_upper_bound(b.var, b.ub)
		end
	end
	nothing
end

"""
    save_bound_if_exists!(reg :: [SavedBound], var :: VariableRef) -> reg

If `var` is fixed or has a lower/upper bound, then a corresponding
SavedBound is pushed into `reg`. The method always return `reg`.
"""
function save_bound_if_exists!(
	reg :: Vector{SavedBound}, var :: VariableRef
) :: Vector{SavedBound}
	was_fixed = is_fixed(var)
	had_lb = has_lower_bound(var)
	had_ub = has_upper_bound(var)
	if was_fixed || had_lb || had_ub
		push!(reg, SavedBound(
			var, was_fixed, (was_fixed ? fix_value(var) : 0.0),
			had_lb, (had_lb ? lower_bound(var) : 0.0),
			had_ub, (had_ub ? upper_bound(var) : 0.0)
		))
	end
	reg
end

"""
    fix_vars!(vars, fix_value = 0.0) :: Vector{SavedBound}

Fix all `vars` to `fix_value` but before save their bound to the returned
vector of `SavedBound`s. The length of the vector returned may be smaller
than the length of `vars` because the bounds are saved only if the variable
actually was fixed or had a bound.
"""
function fix_vars!(vars, fix_value = 0.0)
	saved_bounds = Vector{SavedBound}()
	for var in vars
		save_bound_if_exists!(saved_bounds, var)
	end
	fix.(vars, fix_value; force = true)
	return saved_bounds
end

"""
    relax_vars!(vars, model; check_validity = true) -> tuple of six vectors

The `vars` from `model` are made continuous but if the lower bound is smaller
than zero it becomes zero (the analogue is valid for the upper bound and one).
If `check_validity` is `true` then `!is_valid` variables are skipped.
Four binary vectors (`was_int`, `was_bin`, `had_lb`, `had_ub`) and two
Float64 (`lb_val`, `ub_val`) vectors are returned.

Note: I am not sure anymore why I designed the method to work this way.

See also: [`unrelax_vars!`](@ref), [`relax_all_vars!`](@ref)
"""
function relax_vars!(vars, model; check_validity = true)
	was_int, was_bin, had_lb, had_ub = is_integer.(vars), is_binary.(vars),
		has_lower_bound.(vars), has_upper_bound.(vars)
	lb_val = map(var -> has_lower_bound(var) ? lower_bound(var) : 0.0, vars)
	ub_val = map(var -> has_upper_bound(var) ? upper_bound(var) : 0.0, vars)

	for (i, var) in enumerate(vars)
		check_validity && !is_valid(model, var) && continue
		if was_int[i]
			unset_integer(var)
		elseif was_bin[i]
			unset_binary(var)
			if !had_lb[i] || lb_val[i] < 0.0
				set_lower_bound(var, 0.0)
			end
			if !had_ub[i] || ub_val[i] > 1.0
				set_upper_bound(var, 1.0)
			end
		end
	end

	(was_int, was_bin, had_lb, had_ub, lb_val, ub_val)
end

"""
    relax_vars!(model)

Literally `relax_vars!(all_variables(model), model; check_validity = false)`.
`all_variables` already returns only valid variables.
"""
function relax_all_vars!(model)
	relax_vars!(all_variables(model), model; check_validity = false)
end

"""
    unrelax_vars!(vars, model, original_settings; check_validity = true)

Given the tuple of six vectors returned by `relax_vars!` as
`original_settings`, use it to restore the `vars` from `model` to
their original settings.

See also: [`relax_vars!`](@ref), [`unrelax_all_vars!`](@ref)
"""
function unrelax_vars!(vars, model, original_settings; check_validity = true)
	was_int, was_bin, had_lb, had_ub, lb_val, ub_val = original_settings

	for (i, var) in enumerate(vars)
		check_validity && !is_valid(model, var) && continue
		if was_int[i]
			set_integer(var)
		elseif was_bin[i]
			if !had_lb[i]
				unset_lower_bound(var)
			elseif lb_val[i] != 0.0
				set_lower_bound(var, lb_val[i])
			end
			if !had_ub[i]
				unset_upper_bound(var)
			elseif ub_val[i] != 1.0
				set_upper_bound(var, ub_val[i])
			end
			set_binary(var)
		end
	end

	vars
end

"""
    unrelax_all_vars!(model, original_settings)

If you used `relax_all_vars!` you may use `unrelax_all_vars!` instead of
`unrelax_vars!`.

See also: [`relax_all_vars!`](@ref), [`unrelax_vars!`](@ref)
"""
function unrelax_all_vars!(model, original_settings)
	all_vars = all_variables(model)
	# if variables were deleted since relax_all_vars! the assert will trigger
	@assert length(all_vars) == length(original_settings)
	# the alternatives are: (1) save the list of vars before deletion and pass
	# to unrelax_vars! (checking validity); (2) edit original_settings to remove
	# the references to the settings of deleted vars. we do the first

	# check_validity is false because all_variables only return valid variables
	unrelax_vars!(all_vars, original_settings; check_validity = false)
end

# Not used by the rest of Utilities, just made available to other modules.
include("Args.jl")

end # module
