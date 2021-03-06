#!/bin/bash
# -*- mode: julia -*-
#=
exec julia --project=@. --color=yes --startup-file=no -e "include(popfirst!(ARGS))" "${BASH_SOURCE[0]}" "$@"
=#

using JuMP, Cbc

if length(ARGS) != 1
	print("usage: ./<script_name> <file supported by JuMP.read_from_file>")
end

m = JuMP.read_from_file(ARGS[1])
set_optimizer(m, Cbc.Optimizer)
optimize!(m)
@show JuMP.primal_status(m)
if JuMP.primal_status(m) == MOI.FEASIBLE_POINT
	@show JuMP.objective_value(m)
end

