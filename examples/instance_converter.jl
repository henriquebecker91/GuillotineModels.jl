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

main()

