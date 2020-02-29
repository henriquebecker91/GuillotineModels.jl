using GuillotineModels: CutPattern, to_pretty_str, simplify!

# normalize whitespace
function nws(s :: String) :: String
	strip(replace(s, r"\s+" => ' '))
end

function test_pretty_strs()
	@test to_pretty_str(CutPattern(10, 15, 7)) === "7p10x15"
	@test to_pretty_str(CutPattern(10, 15)) === "P10x15"
	waste = CutPattern(5, 10)
	piece = CutPattern(16, 8, 4)
	h_cut_ex = CutPattern(33, 11, false, [waste, piece])
	@test nws(to_pretty_str(h_cut_ex)) === "P33x11{ P5x10 4p16x8 }"
	h_cut_ex = CutPattern(33, 11, false, [waste, piece])
	@test nws(to_pretty_str(h_cut_ex)) === "P33x11{ P5x10 4p16x8 }"
	patt_and_pieces = [CutPattern(10, 8, 6), h_cut_ex, CutPattern(20, 9, 5)]
	v_cut_ex = CutPattern(33, 28, true, patt_and_pieces)
	v_cut_ex_str = "P33x28[ 6p10x8 P33x11{ P5x10 4p16x8 } 5p20x9 ]"
	@test nws(to_pretty_str(v_cut_ex)) === v_cut_ex_str
end

# TODO: to really test `simplify!` we should create some random solutions
# and check if they keep all the pieces and keep the solution valid, as
# well being the same size or smaller.

test_pretty_strs()

