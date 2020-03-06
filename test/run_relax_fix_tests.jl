import GuillotineModels.Utilities: save_and_relax!, save_and_fix!, restore!
import JuMP
import JuMP: @variable, @constraint, @objective
import GLPK
import MathOptInterface
const MOI = MathOptInterface
import Suppressor: @suppress_out

function test_relax_continuous()
	# Not much to do, just solve an LP and check if the solution keeps being
	# the same after the relax.
	m = JuMP.Model(GLPK.Optimizer)
	@variable(m, 0.5 <= c <= 2)
	@variable(m, -2 <= d <= 0.5)
	# Get the optimal value of the maximization LP.
	@objective(m, Max, c + d)
	JuMP.optimize!(m)
	@test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
	@test JuMP.objective_value(m) ≈ 2.5 rtol=1e-6 atol=1e-6
	# Get the optimal value of the maximization LP (relaxed).
	svc_c = save_and_relax!(c)
	svc_d = save_and_relax!(d)
	JuMP.optimize!(m)
	@test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
	@test JuMP.objective_value(m) ≈ 2.5 rtol=1e-6 atol=1e-6
	# Get the optimal value of the maximization LP (restored).
	restore!.([c, d], [svc_c, svc_d])
	JuMP.optimize!(m)
	@test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
	@test JuMP.objective_value(m) ≈ 2.5 rtol=1e-6 atol=1e-6
end

function test_relax_integer()
	m = JuMP.Model(GLPK.Optimizer)
	@variable(m, 0.5 <= c <= 2, Int)
	@variable(m, -2 <= d <= 0.5, Int)
	# Get the optimal value of the maximization MIP.
	@objective(m, Max, c + d)
	JuMP.optimize!(m)
	@test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
	@test JuMP.objective_value(m) ≈ 2.0 rtol=1e-6 atol=1e-6
	# Get the optimal value of the maximization LP.
	svc_c = save_and_relax!(c)
	svc_d = save_and_relax!(d)
	JuMP.optimize!(m)
	@test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
	@test JuMP.objective_value(m) ≈ 2.5 rtol=1e-6 atol=1e-6
	# Get the optimal value of the minimization LP.
	@objective(m, Min, c + d)
	JuMP.optimize!(m)
	@test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
	@test JuMP.objective_value(m) ≈ -1.5 rtol=1e-6 atol=1e-6
	# Get the optimal value of the restored MIP (but for minimization now).
	restore!.([c, d], [svc_c, svc_d])
	JuMP.optimize!(m)
	@test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
	@test JuMP.objective_value(m) ≈ -1.0 rtol=1e-6 atol=1e-6
end

function test_relax_binary()
	# Given a model with three binary variables, it gives a 2.0 objective value
	# because two variables may only assume value zero.
	m = JuMP.Model(GLPK.Optimizer)
	@variable(m, a, Bin)
	@variable(m, b, Bin)
	@variable(m, 0.5 <= c <= 2, Bin)
	@variable(m, -2 <= d <= 0.5, Bin)
	@objective(m, Max, a + b + c + d)
	@constraint(m, b <= 0.5)
	JuMP.optimize!(m)
	@test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
	@test JuMP.objective_value(m) ≈ 2.0 rtol=1e-6 atol=1e-6

	# If variable `a` is relaxed, the implicit bound at 1 should yet
	# be respected.
	svc_a = save_and_relax!(a)
	JuMP.optimize!(m)
	@test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
	@test JuMP.objective_value(m) ≈ 2.0 rtol=1e-6 atol=1e-6

	# If variable `b` is relaxed, the constraint should keep it at 0.5 anyway.
	svc_b = save_and_relax!(b)
	JuMP.optimize!(m)
	@test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
	@test JuMP.objective_value(m) ≈ 2.5 rtol=1e-6 atol=1e-6

	# If variable `c` is relaxed, the implicit bound at 1 should be respected,
	# even if there is a explicit bound at 2.
	svc_c = save_and_relax!(c)
	JuMP.optimize!(m)
	@test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
	@test JuMP.objective_value(m) ≈ 2.5 rtol=1e-6 atol=1e-6

	# If variable `d` is relaxed, the implicit bound at 0.5 should be respected.
	svc_d = save_and_relax!(d)
	JuMP.optimize!(m)
	@test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
	@test JuMP.objective_value(m) ≈ 3.0 rtol=1e-6 atol=1e-6

	# If the variables are restores to their original state, the same initial
	# objective value should reappear.
	restore!.([a, b, c, d], [svc_a, svc_b, svc_c, svc_d])
	JuMP.optimize!(m)
	@test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
	@test JuMP.objective_value(m) ≈ 2.0 rtol=1e-6 atol=1e-6

	# If the problem minimizes, the implicit bound at c should prevent it
	# becoming zero (instead it stops at 0.5). None should dip below zero.
	@objective(m, Min, a + b + c + d)
	save_and_relax!.([a, b, c, d])
	JuMP.optimize!(m)
	@test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
	@test JuMP.objective_value(m) ≈ 0.5 rtol=1e-6 atol=1e-6
	@test JuMP.value(a) ≈ 0.0 rtol=1e-6 atol=1e-6
	@test JuMP.value(b) ≈ 0.0 rtol=1e-6 atol=1e-6
	@test JuMP.value(c) ≈ 0.5 rtol=1e-6 atol=1e-6
	@test JuMP.value(d) ≈ 0.0 rtol=1e-6 atol=1e-6
end

function test_relax()
	@testset "Relax Binary" begin test_relax_binary() end
	@testset "Relax Integer" begin test_relax_integer() end
	@testset "Relax Continuous" begin test_relax_continuous() end
end

function test_fix()
	m = JuMP.Model(GLPK.Optimizer)
	@variable(m, 10 <= a <= 20)
	@variable(m, b, Bin)
	@variable(m, 30 <= c <= 40, Int)
	@objective(m, Max, a + b + c)
	@constraint(m, b >= 0.1)
	JuMP.optimize!(m)
	@test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
	@test JuMP.objective_value(m) ≈ 61.0 rtol=1e-6 atol=1e-6
	# Fix a continuous variable to a valid value.
	svc_a = save_and_fix!(a, 15.5)
	JuMP.optimize!(m)
	@test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
	@test JuMP.objective_value(m) ≈ 56.5 rtol=1e-6 atol=1e-6
	# Fix the binary variable to a value that would be valid by its binary
	# bounds, but there is a constraint in conflict.
	svc_b = save_and_fix!(b, 0.0)
	JuMP.optimize!(m)
	@test JuMP.primal_status(m) == MOI.NO_SOLUTION
	# Restore the binary variable.
	restore!(b, svc_b)
	JuMP.optimize!(m)
	@test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
	@test JuMP.objective_value(m) ≈ 56.5 rtol=1e-6 atol=1e-6
	# Fix an integer variable to a valid value.
	svc_c = save_and_fix!(c, 35)
	JuMP.optimize!(m)
	@test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
	@test JuMP.objective_value(m) ≈ 51.5 rtol=1e-6 atol=1e-6
	# Fix a integer variable to a value that is not integer.
	svc_c2 = save_and_fix!(c, 35.5)
	@suppress_out JuMP.optimize!(m) # GLPK outputs some warnings if not suppressed
	@test JuMP.primal_status(m) == MOI.NO_SOLUTION
	# Restore the `c` variable from its first save_and_fix! not the last.
	restore!(c, svc_c)
	JuMP.optimize!(m)
	@test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
	@test JuMP.objective_value(m) ≈ 56.5 rtol=1e-6 atol=1e-6
end

test_relax()
test_fix()

