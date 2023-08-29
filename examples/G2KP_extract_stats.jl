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

# The 59 instances presented in table 2.1, page 30,
# from DOI: 10.6092/unibo/amsdottorato/7399 (order was kept).
const THOMOPULOS_THESIS_INSTANCES = vcat(
	# unweighted (37)
	"gcut" .* string.(1:12)
	, String.(split("wang20 2s 3s A1s A2s STS2s STS4s"))
	, ["OF1", "OF2", "W", "CHL1s", "CHL2s"]
	, "A" .* string.(3:5)
	, "CHL" .* string.(5:7)
	, ["CU1", "CU2"]
	, "Hchl" .* split("3s 4s 6s 7s 8s")
	# weighted (22)
	, "cgcut" .* string.(1:3)
	, "okp" .* string.(1:5)
	, String.(split("HH 2 3 A1 A2 STS2 STS4 CHL1 CHL2"))
	, "CW" .* string.(1:3)
	, "Hchl" .* ["2", "9"]
) :: Vector{String}

function main(ARGS = ARGS)
	if iszero(length(ARGS))
		println("USAGE: ./G2KP_extract_stats.jl <Classic_G2KP instance filepaths>")
		exit()
	end
	println("filename;L;W;max_dub;min_dub;qt_piece_types;qt_pieces;total_piece_area;rnl;rnw;rrl;rrw")
	for filepath in ARGS
		filename = basename(filepath)
		data = read_from_file(Val(:Classic_G2KP), filepath)
		L, W = data.L, data.W
		min_d, max_d = minimum(data.d), maximum(data.d)
		qt_piece_types, qt_pieces = length(data.d), sum(data.d)
		total_piece_area = sum(data.l .* data.w .* data.d)
		rnl = rep_number(data.l)
		rnw = rep_number(data.w)
		rrl = rep_ratio(Rational, data.l)
		rrw = rep_ratio(Rational, data.w)
		println("$filename;$L;$W;$max_d;$min_d;$qt_piece_types;$qt_pieces;$total_piece_area;$rnl;$rnw;$rrl;$rrw")
	end
end

#main()

main(ARGS[1] .* THOMOPULOS_THESIS_INSTANCES)

