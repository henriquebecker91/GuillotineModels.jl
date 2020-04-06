#!/bin/bash
# -*- mode: julia -*-
#=
exec julia --project=@. --color=yes --startup-file=no -e "include(popfirst!(ARGS))" "${BASH_SOURCE[0]}" "$@"
=#

import Gurobi, JuMP, GuillotineModels

model = GuillotineModels.CommandLine.SolversArgs.empty_configured_model(
	Val(:Gurobi), Dict{String, Any}([
		"threads" => 1, "no-output" => false, "time-limit" => 3600, "seed" => 0,
		"raw-parameters" => "Pair{String, Any}[]"
	])
)
N, L, W, l, w, p, d = GuillotineModels.InstanceReader.read_from_string(
"""
132 100
20
18 39  702 4
13 49  637 2
27 51 1377 1
65 31 2015 3
45 27 1215 2
69 21 1449 4
21 63 1323 2
54 41 2214 2
41 31 1271 1
37 22  814 1
29 54 1566 1
47 31 1457 2
18 31  558 3
21 17  357 1
19 53 1007 1
13 41  533 4
37 12  444 5
19 29  551 2
67 17 1139 1
49 31 1519 3
"""
)
GuillotineModels.build_model(Val(:PPG2KP), model, d, p, l, w, L, W, Dict{String, Any}(
	"debug" => true, "no-pricing" => true
))

function raw_warm_start(model, nzpe_idxs, nzpe_vals, nzcm_idxs, nzcm_vals)
	@assert length(nzpe_idxs) == length(nzpe_vals)
	@assert length(nzcm_idxs) == length(nzcm_vals)
	pe = model[:picuts]
	cm = model[:cuts_made]
	JuMP.set_start_value.(pe[nzpe_idxs], nzpe_vals)
	JuMP.set_start_value.(cm[nzcm_idxs], nzcm_vals)
end

raw_warm_start(model,
#=
	[704, 740, 813, 849, 1158, 1251, 1294, 1409, 1412, 1930, 3072, 3832],
	[1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1],
	[39, 2207, 3998, 27744, 34512, 34721, 35706, 36807, 37412, 38732, 42992],
	[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
=#
	[704, 740, 813, 849, 1158, 1251, 1294, 1409, 1412, 1930, 3072, 3832],
	[1, 1, 1, 1, 0, 1, 1, 3, 0, 1, 0, 1],
	[39, 2207, 3998, 27744, 34512, 34721, 35706, 36807, 37412, 38732, 42992],
	[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
)

JuMP.optimize!(model)

