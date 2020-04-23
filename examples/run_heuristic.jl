#!/bin/bash
# -*- mode: julia -*-
#=
exec julia --project=@. --color=yes --startup-file=no -e "include(popfirst!(ARGS))" "${BASH_SOURCE[0]}" "$@"
=#

import GuillotineModels.InstanceReader.read_from_file
using GuillotineModels.PPG2KP.Heuristic
using TimerOutputs
using RandomNumbers.Xorshifts # for Xoroshiro128Plus (for heuristic rng)

function main(args)
	timer = TimerOutput()
	# Parsing.
	qt_args = length(args)
	if qt_args < 1 || qt_args > 3
		println("usage: ./run_heuristic instance_filename [rng_seed] [max_iter]")
		exit()
	end
	filename = args[1]
	seed = (qt_args > 1 ? args[2] : 1)
	max_iter = (qt_args > 2 ? args[3] : 1000000)

	# Initializing.
	N, L, W, l, w, p, d = read_from_file(filename)

	GC.enable(false)
	# Running.
	rng = Xoroshiro128Plus(seed)
	bkv2, sel2, pat2 = @timeit timer "mock_fast" begin
		bkv, sel, pat = fast_iterated_greedy(
			d, p, l, w, L, W, rng, 1
		)
		GC.gc()
		bkv, sel, pat
	end
	rng = Xoroshiro128Plus(seed)
	bkv1, sel1, pat1 = @timeit timer "mock_normal" begin
		bkv, sel, pat = iterated_greedy(
			d, p, l, w, L, W, rng, 1
		)
		GC.gc()
		bkv, sel, pat
	end
	if (bkv1, sel1, pat1) == (bkv2, sel2, pat2)
		println("The two mocks give the same values.")
		println("bkv = $(bkv1)")
		println("sel = $(sel1)")
		println("pat = $(pat1)")
	else
		println("The two mocks give distinct values. They follow:")
		println("iterated_greedy return:")
		@show bkv1
		@show sel1
		@show pat1
		println("fast_iterated_greedy return:")
		@show bkv2
		@show sel2
		@show pat2
	end

	# TODO: disable GC and the re-enable it and call it for each call
	rng = Xoroshiro128Plus(seed)
	bkv1, sel1, pat1 = @timeit timer "normal" begin
		bkv, sel, pat = iterated_greedy(
			d, p, l, w, L, W, rng, max_iter
		)
		GC.gc()
		bkv, sel, pat
	end
	rng = Xoroshiro128Plus(seed)
	bkv2, sel2, pat2 = @timeit timer "fast" begin
		bkv, sel, pat = fast_iterated_greedy(
			d, p, l, w, L, W, rng, max_iter
		)
		GC.gc()
		bkv, sel, pat
	end
	if (bkv1, sel1, pat1) == (bkv2, sel2, pat2)
		println("The two runs give the same values. They follow:")
		println("bkv = $(bkv1)")
		println("sel = $(sel1)")
		println("pat = $(pat1)")
	else
		println("The two runs give distinct values. They follow:")
		println("iterated_greedy return:")
		@show bkv1
		@show sel1
		@show pat1
		println("fast_iterated_greedy return:")
		@show bkv2
		@show sel2
		@show pat2
	end
	TimerOutputs.print_timer(timer)
	println()
end

main(ARGS)

