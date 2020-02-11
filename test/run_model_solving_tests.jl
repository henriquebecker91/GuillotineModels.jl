const TINY_HANDMADE_INSTANCES = Tuple{String, Int64, String}[
	("Plate 1x1 and 1 piece 1x1 with profit 1.", 1, """
		1 1
		1
		1 1 1 1
	"""),
	("Pieces do not fit plate.", 0, """
		10 10
		3
		15 15 1 1
		15 10 1 1
		10 15 1 1
	"""),
	("One piece fit the plate others do not.", 10, """
		10 10
		4
		10 10 10 1
		15 15 1 1
		15 10 1 1
		10 15 1 1
	""")
]

const EASY_LITERATURE_INSTANCES = Tuple{String, Int64, String}[
  ("gcut1", 48368, """
		250 250
		10
		167 184 30728 1
		114 118 13452 1
		167 152 25384 1
		83 140 11620 1
		70 86 6020 1
		143 166 23738 1
		120 160 19200 1
		66 148 9768 1
		87 141 12267 1
		69 165 11385 1
	""")
]

using Test
import JuMP
import MathOptInterface
const MOI = MathOptInterface
import GLPK # to load glue code in GuillotineModels
import Cbc
using GuillotineModels: build_model
using GuillotineModels.InstanceReader: read_from_string
using GuillotineModels.Utilities.Args: accepted_arg_list,
	create_normalized_arg_subset
using GuillotineModels.CommandLine.SolversArgs: empty_configured_model

function test_obj_val(
	model :: Symbol, solver :: Symbol, instance :: String, obj_val :: Int64
)
	N, L, W, l, w, p, d = read_from_string(instance)
	solver_type = Val(solver)
	solver_pp = create_normalized_arg_subset(
		Dict{String, Any}("no-output" => true),
		accepted_arg_list(solver_type)
	)
	m = empty_configured_model(solver_type, solver_pp)
	build_model(Val(model), m, d, p, l, w, L, W, Dict{String, Any}())
	JuMP.optimize!(m)
	@test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
	@test JuMP.objective_value(m) == obj_val
end

function test_obj_val_of_all_combinations(
	models :: Vector{Symbol},
	solvers :: Vector{Symbol},
	instances :: Vector{Tuple{String, Int64, String}}
)
	@testset "Model: $model" for model in models
		@testset "Solver: $solver" for solver in solvers
			for (description, obj_val, instance) in instances
				@testset "Instance: $description" begin
					test_obj_val(model, solver, instance, obj_val)
				end
			end
		end
	end
end

test_obj_val_of_all_combinations(
	[:PPG2KP, :Flow],
	[:GLPK],
	[TINY_HANDMADE_INSTANCES; EASY_LITERATURE_INSTANCES]
)

