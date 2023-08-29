#!/bin/bash
# -*- mode: julia -*-
#=
exec julia --compile=min -O0 --project=@. --color=yes --startup-file=no -e "include(popfirst!(ARGS))" "${BASH_SOURCE[0]}" "$@"
=#

using GuillotineModels.Data

function rep_number(type, xs)
	seen = Set{eltype(xs)}()
	qt_repetitions = zero(type)
	for x in xs
		if x in seen
			qt_repetitions += one(type)
		else
			push!(seen, x)
		end
	end
	return qt_repetitions
end

rep_number(xs) = rep_number(typeof(length(xs)), xs)

function rep_ratio(type, xs)
	max_repetitions = convert(type, length(xs)) - one(type)
	return rep_number(type, xs) / max_repetitions
end

rep_ratio(xs) = rep_ratio(typeof(length(xs)), xs)

function main(ARGS = ARGS)
	if iszero(length(ARGS))
		println("USAGE: ./CPG_SSSCSP_extract_stats.jl <CPG_SSSCSP instance filepaths>")
		exit()
	end
	println("filename;L;W;max_dub;min_dub;qt_piece_types;qt_pieces;total_piece_area;rnl;rnw;rrl;rrw")
	for filepath in ARGS
		filename = basename(filepath)
		data = only(read_from_file(Val(:CPG_SSSCSP), filepath))
		L, W = data.L, data.W
		min_dub, max_dub = minimum(data.d), maximum(data.d)
		qt_piece_types, qt_pieces = length(data.d), sum(data.d)
		total_piece_area = sum(data.l .* data.w .* data.d)
		rnl = rep_number(data.l)
		rnw = rep_number(data.w)
		rrl = rep_ratio(Rational, data.l)
		rrw = rep_ratio(Rational, data.w)
		println("$filename;$L;$W;$max_dub;$min_dub;$qt_piece_types;$qt_pieces;$total_piece_area;$rnl;$rnw;$rrl;$rrw")
	end
end

main()

