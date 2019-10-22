module ModelUtils

using JuMP

export num_all_constraints, reduced_cost, delete_vars_by_pricing!
export relax_vars!, relax_all_vars!, unrelax_vars!, unrelax_all_vars!
export SavedBound, save_bound_if_exists!, restore_bound!
export which_vars_to_keep, fix_vars!

# JuMP has no method for getting all constraints, you need to get the
# types of constraints used in the model and then query the number for
# specific type.
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

# NOTE: by design, this method do not change the vars vector itself, but
# instead just calls delete over part of its contents (what mark such content
# as invalid)
# NOTE: this method was abandoned because deleting 4
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

# Execute Furini's 2016 Final Pricing. The returned value is a boolean list,
# with each index corresponding to a variable in all_variables(model),
# and with a true value if the variable should be kept, and false if it should
# be removed. The model parameter is expected to be a solved continuous
# relaxation of the problem (so the objective_value and the reduced_cost may
# be extracted).
function which_vars_to_keep(model, lb :: Float64)
  obj = objective_value(model)
  rcs = reduced_cost.(all_variables(model))
  # We keep the ones equal to lb (not strictly necessary) to guarantee the
  # warm start will work (the variables needed for it should be there).
  return map(rc -> obj - rc >= lb, rcs)
end

struct SavedBound
  var :: VariableRef
  was_fixed :: Bool
  fix_value :: Float64
  had_lb :: Bool
  lb :: Float64
  had_ub :: Bool
  ub :: Float64
end

function restore_bound!(b :: SavedBound) :: Nothing
  if b.was_fixed && !(is_fixed(b.var) && fix_value(b.var) == b.fix_value)
    fix(b.var, b.fix_value; force = true)
  else
    if b.has_lb && !(has_lower_bound(b.var) && lower_bound(b.var) == b.lb)
      set_lower_bound(b.var, b.lb; force = true)
    end
    if b.has_ub && !(has_upper_bound(b.var) && upper_bound(b.var) == b.ub)
      set_upper_bound(b.var, b.ub; force = true)
    end
  end
  nothing
end

function save_bound_if_exists!(
  reg :: Vector{SavedBound}, var :: VariableRef
) :: Vector{SavedBound}
  was_fixed = is_fixed(var)
  had_lb = has_lower_bound(var)
  had_ub = has_upper_bound(var)
  if was_fixed || had_lb || had_ub
    push!(reg, SavedBound(
      var, was_fixed, was_fixed ? fix_value(var) : 0.0,
      had_lb, had_lb ? lower_bound(var) : 0.0,
      had_ub, had_ub ? upper_bound(var) : 0.0
    ))
  end
  reg
end

function fix_vars!(vars, fix_value = 0.0)
  saved_bounds = Vector{SavedBound}()
  for var in vars
    save_bound_if_exists!(saved_bounds, var)
  end
  fix.(vars, fix_value; force = true)
  return saved_bounds
end

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

function relax_all_vars!(model)
  # all_variables already just return valid variables
  relax_vars!(all_variables(model), model; check_validity = false)
end

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

end # module
