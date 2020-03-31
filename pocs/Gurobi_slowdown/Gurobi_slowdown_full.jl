#!/bin/bash
# -*- mode: julia -*-
#=
exec julia --project=@. --color=yes --startup-file=no -e "include(popfirst!(ARGS))" "${BASH_SOURCE[0]}" "$@"
=#

using JuMP
import Gurobi

function alternating()
    m = direct_model(Gurobi.Optimizer())
    @variable(m, x[1:100000] >= 0, Int)
    for e in x
        unset_integer(e)
        is_integer(e)
        lower_bound(e)
    end
    return
end

function batch()
    m = direct_model(Gurobi.Optimizer())
    @variable(m, x[1:100000] >= 0, Int)
    unset_integer.(x)
    is_integer.(x)
    lower_bound.(x)
    return
end

struct SavedVarConf
  was_bin :: Bool
  was_int :: Bool
  was_fixed :: Bool
  fix_value :: Float64
  had_lb :: Bool
  lb :: Float64
  had_ub :: Bool
  ub :: Float64
end

function SavedVarConf(var :: VariableRef) :: SavedVarConf
  was_bin = is_binary(var)
  was_int = is_integer(var)
  was_fixed = is_fixed(var)
  had_lb = has_lower_bound(var)
  had_ub = has_upper_bound(var)
  return SavedVarConf(
    was_bin, was_int,
    was_fixed, (was_fixed ? fix_value(var) : 0.0),
    had_lb, (had_lb ? lower_bound(var) : 0.0),
    had_ub, (had_ub ? upper_bound(var) : 0.0)
  )
end

function relax!(var :: VariableRef) :: VariableRef
  if is_integer(var)
    unset_integer(var)
  elseif is_binary(var)
    unset_binary(var)
    if !has_lower_bound(var) || lower_bound(var) < 0.0
      set_lower_bound(var, 0.0)
    end
    if !has_upper_bound(var) || upper_bound(var) < 0.0
      set_upper_bound(var, 1.0)
    end
  end

  return var
end

function save_and_relax!(var :: VariableRef)
  svc = SavedVarConf(var)
  if svc.was_int
    unset_integer(var)
  elseif svc.was_bin
    unset_binary(var)
    if !svc.had_lb || svc.lb < 0.0
      set_lower_bound(var, 0.0)
    end
    if !svc.had_ub || svc.ub > 1.0
      set_upper_bound(var, 1.0)
    end
  end

  return svc
end

function alternating2()
    m = direct_model(Gurobi.Optimizer())
    @variable(m, x[1:100000] >= 0, Int)
    save_and_relax!.(x)
    return
end

function batch2()
    m = direct_model(Gurobi.Optimizer())
    @variable(m, x[1:100000] >= 0, Int)
    SavedVarConf.(x)
    relax!.(x)
    return
end

function test()
    println("-------------------- batch --------------------------")
    @time batch()
    println("-------------------- alternating --------------------")
    @time alternating()
    #=
    println("-------------------- batch2 -------------------------")
    @time batch2()
    println("-------------------- alternating2 -------------------")
    @time alternating2()
    =#
end

test()
#println("-------------------- second run ---------------------")
#test()

