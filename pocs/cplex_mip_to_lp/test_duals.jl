import CPLEX, Gurobi, GLPK
using JuMP

function test_duals(solver_sym, relaxed_mip)
	solver_lib = getfield(Main, solver_sym)
	m = direct_model(solver_lib.Optimizer())
	set_silent(m)
	if relaxed_mip
		@variable(m, x[1:10] >= 0, Int)
	else
		@variable(m, x[1:10] >= 0)
	end
	@objective(m, Min, sum(x))
	relaxed_mip && unset_integer.(x)
	optimize!(m)
	@show solver_sym
	@show relaxed_mip
	@show has_duals(m)
end

test_duals(:GLPK, false)
test_duals(:GLPK, true)
test_duals(:CPLEX, false)
test_duals(:CPLEX, true)
test_duals(:Gurobi, false)
test_duals(:Gurobi, true)
