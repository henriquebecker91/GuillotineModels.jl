const TINY_HANDMADE_INSTANCES = Tuple{String, Int64, String}[
	("Plate 1x1 and 1 piece 1x1 with profit 1.", 1, """
		1 1
		1
		1 1 1 1
	""")
	, ("Pieces do not fit plate.", 0, """
		10 10
		3
		15 15 1 1
		15 10 1 1
		10 15 1 1
	""")
	, ("One piece fit the plate others do not.", 10, """
		10 10
		4
		10 10 10 1
		15 15 1 1
		15 10 1 1
		10 15 1 1
	""")
	, ("Two halfs instead slightly larger very efficient piece.", 20, """
		10 20
		2
		10 10 10 2
		10 11 19 1
	""")
	, ("The whole plate could be used if the plates were rotated.", 20, """
		16 20
		2
		10 8 10 2
		10 16 19 1
	""")
	, ("Plate 100x100 and 100 pieces 10x10 with profit 7.", 700, """
		100 100
		1
		10 10 7 100
	""")
	, ("Plate 100x100 and 50 pieces 10x10 with profit 7.", 350, """
		100 100
		1
		10 10 7 50
	""")
	, ("Plate 50x200 and 200 pieces 10x10 with profit 7.", 700, """
		100 100
		1
		10 10 7 200
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
	, ("cgcut1", 244, """
		15 10
		7
		8 4 66 2
		3 7 35 1
		8 2 24 3
		3 4 17 5
		3 3 11 2
		3 2  8 2
		2 1  2 1
	""")
	, ("CHL5", 390, """
		20 20
		10
		14  2 28 1
		1   5  5 1
		20  4 80 2
		12  3 36 3
		11  8 88 2
		11  6 66 2
		7   9 63 1
		17  5 85 1
		7  14 98 2
		1   7  7 3
	""")
]

import JuMP
import MathOptInterface
const MOI = MathOptInterface
import GLPK # to load glue code in GuillotineModels
import Cbc
using GuillotineModels: build_model, get_cut_pattern, to_pretty_str, simplify!
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
	bmr = build_model(Val(model), m, d, p, l, w, L, W, Dict{String, Any}())
	#JuMP.write_to_file(m, "Cbc_false_infeasible.lp")
	#JuMP.write_to_file(m, "Cbc_false_infeasible.mof.json")
	#JuMP.write_to_file(m, "Cbc_false_infeasible.mps")
	JuMP.optimize!(m)
	@test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
	@test JuMP.objective_value(m) ≈ obj_val rtol=1e-6 atol=1e-6
	#=
	if model === :PPG2KP
		println(instance)
		solution = get_cut_pattern(Val(model), m, eltype(d), eltype(l), bmr)
		@show solution
		println("to_pretty_str(solution) = $(to_pretty_str(solution))")
		println("to_pretty_str(simplify!(solution)) = $(to_pretty_str(simplify!(solution)))")
	end
	=#
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
					#println(description)
					test_obj_val(model, solver, instance, obj_val)
				end
			end
		end
	end
end

test_obj_val_of_all_combinations(
	[:PPG2KP, :Flow],
	[:GLPK, #=:Cbc=#], # Cbc is kinda bugged
	[TINY_HANDMADE_INSTANCES; EASY_LITERATURE_INSTANCES]
)

