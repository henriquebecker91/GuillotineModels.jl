#!/bin/bash
# -*- mode: julia -*-
#=
exec julia --compile=min -O0 --project=@. --color=yes --startup-file=no -e "include(popfirst!(ARGS))" "${BASH_SOURCE[0]}" "$@"
=#

using GuillotineModels.Data: read_from_file, write_to_file, is_collection

function number_path(s, i)
	if match(r"\.", s) !== nothing
		replace(s, '.' => "_$i."; count = 1)
	else
		s * "_$i"
	end
end

function main(ARGS = ARGS)
	if length(ARGS) != 4
		println("USAGE: ./instance_converter.jl <input format> <input file> <output format> <output file>")
		println("If the input format has multiple instances, and the output format just support a single instance, then an output file will be created for each instance and '_X' (where X is the instance number) will be appended before the first dot in the output filename (or just before the end, if necessary).")
		exit()
	end
	in_format, in_file, out_format, out_file = ARGS
	in_format_val = Val(Symbol(in_format))
	instance = read_from_file(in_format_val, in_file)
	out_format_val = Val(Symbol(out_format))
	if is_collection(in_format_val) && !is_collection(out_format_val)
		for (index, inst) in enumerate(instance)
			write_to_file(out_format_val, inst, number_path(out_file, index))
		end
	else
		write_to_file(out_format_val, instance, out_file)
	end
end

function furini2016_convert(input_path, output_path)
	THOMOPULOS_THESIS_INSTANCES = vcat(
		# unweighted
		"gcut" .* string.(1:12)
		, String.(split("wang20 2s 3s A1s A2s STS2s STS4s"))
		, ["OF1", "OF2", "W", "CHL1s", "CHL2s"]
		, "A" .* string.(3:5)
		, "CHL" .* string.(5:7)
		, ["CU1", "CU2"]
		, "Hchl" .* split("3s 4s 6s 7s 8s")
		# weighted
		, "cgcut" .* string.(1:3)
		, "okp" .* string.(1:5)
		, String.(split("HH 2 3 A1 A2 STS2 STS4 CHL1 CHL2"))
		, "CW" .* string.(1:3)
		, "Hchl" .* ["2", "9"]
	) :: Vector{String}
	for inst_name in THOMOPULOS_THESIS_INSTANCES
		main([
			"Classic_G2KP",
			joinpath(input_path, inst_name),
			"Simple_CPG_SLOPP",
			joinpath(output_path, inst_name)
		])
	end
end

main()

# Comment the `main()` call and uncomment the below to have the script
# making all 59 Furini2016 instances (from Classic_G2KP format) into the
# Simple_CPG_SLOPP format.

#furini2016_convert(ARGS...)

