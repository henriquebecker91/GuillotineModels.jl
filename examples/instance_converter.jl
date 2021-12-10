#!/bin/bash
# -*- mode: julia -*-
#=
exec julia --compile=min -O0 --project=@. --color=yes --startup-file=no -e "include(popfirst!(ARGS))" "${BASH_SOURCE[0]}" "$@"
=#

using GuillotineModels.Data
# write_to_file is not exported to avoid clashes with MathOptInterface.
# We also need import instead of using because we gonna extend it.
import GuillotineModels.Data: write_to_file

const SLOPP_HEADER = """
***2D Rectangular Problem***
***Instances for the Single Large Object Placement Problem (SLOPP)***
Input parameter file: instance_not_generated_but_converted_from_another_format
***************************************************************************************************************
Total number of instances 
LargeObject.Length	LargeObject.Width
Number of different item types (i)
Item[i].Length	Item[i].Width	Item[i].LowerBoundDemand	Item[i].UpperBoundDemand	Item[i].Value
***************************************************************************************************************"""

# Extends the GuillotineModels.Data.write_to_file to support writing
# G2KP objects in SLOPP format. Note that below we call
# `write_to_file(format, instance, *path_str*)` but we do not need to
# define the method that takes filenames instead of IO objects because
# the fallback method in GuillotineModels.Data already does this for us.
function write_to_file(
	format :: CPG_SLOPP{D, S, P},
	instance :: G2KP{D, S, P},
	io :: IO
) where {D, S, P}
	L, W, l, w, d, p = instance.L, instance.W, instance.l, instance.w,
		instance.d, instance.p
	write(io, "$SLOPP_HEADER\n")
	write(io, "1\n")
	write(io, "$L $W\n")
	write(io, "$(length(p))\n")
	for (pl, pw, pd, pp) in zip(l, w, d, p)
		write(io, "$pl $pw 0 $pd $pp\n")
	end
	flush(io)
	return io
end

# We also define a default version in which the user only gives the
# format name (wrapped in a Val) and the type of Demand, Size, and
# Profit is guessed from the input instance.
function write_to_file(
	_ :: Val{:CPG_SLOPP},
	instance :: G2KP{D, S, P},
	io :: IO
) where {D, S, P}
	write_to_file(CPG_SLOPP{D, S, P}(), instance, io)
end

function main(ARGS = ARGS)
	if length(ARGS) != 2
		println("USAGE: ./instance_converter.jl <Classic_G2KP input filepath> <CPG_SLOPP output filepath>")
		exit()
	end
	in_file, out_file = ARGS
	in_format_val = Val(:Classic_G2KP)
	instance = read_from_file(in_format_val, in_file)
	out_format_val = Val(:CPG_SLOPP)
	@assert !is_collection(in_format_val)
	write_to_file(out_format_val, instance, out_file)
end

main()

