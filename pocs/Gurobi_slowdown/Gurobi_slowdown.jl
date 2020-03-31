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
        #is_integer(e)
        #has_upper_bound(e)
        lower_bound(e)
    end
    return
end

function batch()
    m = direct_model(Gurobi.Optimizer())
    @variable(m, x[1:100000] >= 0, Int)
    unset_integer.(x)
    is_integer.(x)
    #is_integer.(x)
    #has_upper_bound.(x)
    lower_bound.(x)
    return
end

function test()
    println("-------------------- batch --------------------------")
    @time batch()
    println("-------------------- alternating --------------------")
    @time alternating()
end

test()
println("\n-------------------- second run ---------------------\n")
test()

# The output of this code in my notebook (Intel(R) Core(TM) i7-7700HQ
# CPU @ 2.80GHz) before PR https://github.com/JuliaOpt/Gurobi.jl/pull/303
# was:
#
# -------------------- batch --------------------------
# Academic license - for non-commercial use only
#   0.758784 seconds (2.34 M allocations: 99.862 MiB, 16.04% gc time)
# -------------------- alternating --------------------
# Academic license - for non-commercial use only
#  24.357266 seconds (1.20 M allocations: 47.019 MiB, 0.08% gc time)
#
# -------------------- second run ---------------------
#
# -------------------- batch --------------------------
# Academic license - for non-commercial use only
#   0.210978 seconds (1.31 M allocations: 50.068 MiB, 8.00% gc time)
# -------------------- alternating --------------------
# Academic license - for non-commercial use only
#  23.389932 seconds (1.20 M allocations: 47.019 MiB, 0.04% gc time)


