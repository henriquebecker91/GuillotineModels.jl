#!/bin/bash
# -*- mode: julia -*-
#=
exec julia --compile=min -O0 --project=@. --color=yes --startup-file=no -e "include(popfirst!(ARGS))" "${BASH_SOURCE[0]}" "$@"
=#

using GuillotineModels.Data

function main(ARGS = ARGS)
	if iszero(length(ARGS))
		println("USAGE: ./CPG_MHLOPPW_extract_stats.jl <CPG_MHLOPPW instance filepaths>")
		exit()
	end
	println("filename;L;W;qt_original_plates;max_dub;min_dub;qt_piece_types;qt_pieces;total_piece_area")
	for filepath in ARGS
		filename = basename(filepath)
		data = only(read_from_file(Val(:CPG_MHLOPPW), filepath))
		L, W = only(data.L), only(data.W)
		qt_original_plates = only(data.available)
		min_dub, max_dub = minimum(data.dub), maximum(data.dub)
		qt_piece_types, qt_pieces = length(data.dub), sum(data.dub)
		total_piece_area = sum(data.l .* data.w .* data.dub)
		println("$filename;$L;$W;$qt_original_plates;$max_dub;$min_dub;$qt_piece_types;$qt_pieces;$total_piece_area")
	end
end

main()

