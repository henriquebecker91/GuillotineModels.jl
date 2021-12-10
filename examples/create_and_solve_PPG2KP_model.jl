#!/bin/bash
# -*- mode: julia -*-
#=
exec julia --project=@. --color=yes --startup-file=no -e "include(popfirst!(ARGS))" "${BASH_SOURCE[0]}" "$@"
=#

import JuMP
import GLPK
import GuillotineModels
const GM = GuillotineModels
import MathOptInterface
const MOI = MathOptInterface

function create_and_solve_PPG2KP_model(g2kp_inst_str)
	model = JuMP.Model(GLPK.Optimizer)

	instance = GM.Data.read_from_string(Val(:Classic_G2KP), g2kp_inst_str)
	# bsr stands for build stop reason. mbp stands for model byproduct.
	bsr, mbp = GM.build_model(
		Val(:G2KP), Val(:PPG2KP), instance, model,
		Dict{String,Any}("faithful2furini2016" => true)
	)
	# As nor pricing, nor warm-starting, is enabled. There is no way
	# build_model can find the optimal during the building process.
	@assert bsr == GM.BUILT_MODEL

	JuMP.optimize!(model)

	@show JuMP.termination_status(model)
	if JuMP.primal_status(model) == MOI.FEASIBLE_POINT
		@show JuMP.objective_value(model)
		cp = GM.get_cut_pattern(Val(:G2KP), Val(:PPG2KP), model, mbp)
		println("Human and computer-friendly representation")
		println(GM.to_pretty_str(cp))
		println("Tikz code to plot the solution")
		println(GM.to_tikz_picture(cp))
		GM.simplify!(cp)
		println("Human and computer-friendly representation (simplified solution)")
		println(GM.to_pretty_str(cp))
		println("Tikz code to plot the simplified solution")
		println(GM.to_tikz_picture(cp))
	end
end

# cgcut1 instance
const instance_str = """
15 10
7
8 4 66 2
3 7 35 1
8 2 24 3
3 4 17 5
3 3 11 2
3 2  8 2
2 1  2 1"""

create_and_solve_PPG2KP_model(instance_str)

