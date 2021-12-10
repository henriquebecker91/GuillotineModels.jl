#!/bin/bash
# -*- mode: julia -*-
#=
exec julia --project=@. --color=yes --startup-file=no -e "include(popfirst!(ARGS))" "${BASH_SOURCE[0]}" "$@"
=#

import GuillotineModels

import GuillotineModels.Data: read_from_string, Classic_G2KP
import GuillotineModels.PPG2KP.Enumeration: gen_cuts
#const E = GuillotineModels.PPG2KP.Enumeration
import GuillotineModels.Utilities: SortedLinkedLW

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

let # Makes the scope local avoiding performance loss of global scope.
# Reads the instance and use a different integer type for
# Demand, Size, and Profit/Area.
instance = read_from_string(Classic_G2KP{Int8, Int16, Int32}(), instance_str)
# Creates auxiliar data structure necessary for gen_cuts.
sllw = SortedLinkedLW(Int8, instance.l, instance.w)
# The gen_cuts
byproduct = gen_cuts(
	Int32, instance.d, sllw, instance.L, instance.W; faithful2furini2016 = true
)
println("Each tuple in pli_lwb has length, width, and bound (i.e., max copies) of an enumerated plate.")
@show byproduct.pli_lwb
println("Each tuple in cuts has the index of the parent plate and the indexes of both child plates (an index of zero represents the child is waste).")
@show byproduct.cuts
println("Below we combine the two and give all cuts in terms of size of the parent plate and size of the childs.")
for cut in byproduct.cuts
	pp, fc, sc = cut # parent plate, first child, second child
	print(byproduct.pli_lwb[pp][1], "x", byproduct.pli_lwb[pp][2])
	print(" -> ")
	if !iszero(fc)
		print(byproduct.pli_lwb[fc][1], "x", byproduct.pli_lwb[fc][2])
	end
	if !iszero(sc)
		!iszero(fc) && print(" ")
		print(byproduct.pli_lwb[sc][1], "x", byproduct.pli_lwb[sc][2])
	end
	println()
end

end # let

