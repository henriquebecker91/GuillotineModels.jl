#!/usr/bin/env julia --project=@.

using JuMP, Cbc
#using JuMP, GLPK

function run()
	m = Model(Cbc.Optimizer)
	#m = Model(GLPK.Optimizer)

	empty_terms = AffExpr()
	@constraint(m, empty_terms <= 2)
	optimize!(m)
	@show JuMP.primal_status(m)
	if JuMP.primal_status(m) == MOI.FEASIBLE_POINT
		@show JuMP.objective_value(m)
	end
end

@time run()
@time run()

