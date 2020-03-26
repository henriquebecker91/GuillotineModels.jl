#!/bin/bash
# -*- mode: julia -*-
#=
exec julia --project=@. --color=yes --startup-file=no -e "include(popfirst!(ARGS))" "${BASH_SOURCE[0]}" "$@"
=#

using JuMP, Cbc, GLPK

function run(optimizer)
	m = Model(optimizer)

	empty_terms = AffExpr()
	@constraint(m, empty_terms <= 2)
	optimize!(m)
	@show JuMP.primal_status(m)
	if JuMP.primal_status(m) == MOI.FEASIBLE_POINT
		@show JuMP.objective_value(m)
	end
end

println("\n===================== GLPK First Time ====================\n")
@time run(GLPK.Optimizer)
println("\n==================== GLPK Second Time ====================\n")
@time run(GLPK.Optimizer)
println("\n====================== Cbc First Time ====================\n")
@time run(Cbc.Optimizer)
println("\n===================== Cbc Second Time ====================\n")
@time run(Cbc.Optimizer)

