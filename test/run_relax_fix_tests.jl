import GuillotineModels.Utilities: relax!, restore!
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
	svc_c = relax!(c)
	svc_d = relax!(d)
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
	svc_c = relax!(c)
	svc_d = relax!(d)
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
	svc_a = relax!(a)
	JuMP.optimize!(m)
	@test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
	@test JuMP.objective_value(m) ≈ 2.0 rtol=1e-6 atol=1e-6

	# If variable `b` is relaxed, the constraint should keep it at 0.5 anyway.
	svc_b = relax!(b)
	JuMP.optimize!(m)
	@test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
	@test JuMP.objective_value(m) ≈ 2.5 rtol=1e-6 atol=1e-6

	# If variable `c` is relaxed, the implicit bound at 1 should be respected,
	# even if there is a explicit bound at 2.
	svc_c = relax!(c)
	JuMP.optimize!(m)
	@test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
	@test JuMP.objective_value(m) ≈ 2.5 rtol=1e-6 atol=1e-6

	# If variable `d` is relaxed, the implicit bound at 0.5 should be respected.
	svc_d = relax!(d)
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
	relax!.([a, b, c, d])
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

test_relax()

