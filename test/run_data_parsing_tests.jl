using GuillotineModels.Data

const SLOPP_HEADER = """
***2D Rectangular Problem***
***Instances for the Single Large Object Placement Problem (SLOPP)***
Input parameter file: SLOPP_parameters.txt
***************************************************************************************************************
Total number of instances
LargeObject.Length      LargeObject.Width
Number of different item types (i)
Item[i].Length  Item[i].Width   Item[i].LowerBoundDemand        Item[i].UpperBoundDemand        Item[i].Value
***************************************************************************************************************"""

const SLOPP_BODY_1 = """
50      60
28
21      13      4       4       264
6       45      1       2       191
5       60      2       4       264
30      39      1       1       795
14      42      5       5       311
29      17      2       3       314
10      9       4       5       73
23      31      2       3       567
16      22      4       5       181
14      26      1       1       182
24      24      2       5       476
12      56      2       2       504
12      49      3       5       337
12      47      1       3       413
18      32      4       4       520
26      42      1       2       1064
14      16      2       3       163
8       38      4       5       219
6       58      3       4       293
21      40      2       3       791
13      52      2       5       479
8       17      1       1       73
29      56      1       5       1449
18      9       1       5       101
27      27      1       1       596
6       50      1       2       255
38      21      3       4       571
27      35      2       3       697"""

SLOPP_BODY_2 = """
83      67
10
48      28      4       4       1017
68      44      1       1       2458
14      58      4       5       519
5       43      1       1       205
63      8       1       2       349
49      56      1       5       1447
6       22      3       3       95
12      30      4       4       239
50      42      4       5       1763
17      34      4       4       325 """

const SLOPP_INST_1 = SLOPP{Int,Int,Int}(
	50, 60,
	[21, 6, 5, 30, 14, 29, 10, 23, 16, 14, 24, 12, 12, 12, 18, 26, 14, 8, 6, 21, 13, 8, 29, 18, 27, 6, 38, 27],
	[13, 45, 60, 39, 42, 17, 9, 31, 22, 26, 24, 56, 49, 47, 32, 42, 16, 38, 58, 40, 52, 17, 56, 9, 27, 50, 21, 35],
	[4, 1, 2, 1, 5, 2, 4, 2, 4, 1, 2, 2, 3, 1, 4, 1, 2, 4, 3, 2, 2, 1, 1, 1, 1, 1, 3, 2],
	[4, 2, 4, 1, 5, 3, 5, 3, 5, 1, 5, 2, 5, 3, 4, 2, 3, 5, 4, 3, 5, 1, 5, 5, 1, 2, 4, 3],
	[264, 191, 264, 795, 311, 314, 73, 567, 181, 182, 476, 504, 337, 413, 520, 1064, 163, 219, 293, 791, 479, 73, 1449, 101, 596, 255, 571, 697]
)

const SLOPP_INST_2 = SLOPP{Int64,Int64,Int64}(
	83, 67,
	[48, 68, 14, 5, 63, 49, 6, 12, 50, 17],
	[28, 44, 58, 43, 8, 56, 22, 30, 42, 34],
	[4, 1, 4, 1, 1, 1, 3, 4, 4, 4],
	[4, 1, 5, 1, 2, 5, 3, 4, 5, 4],
	[1017, 2458, 519, 205, 349, 1447, 95, 239, 1763, 325]
)

const SLOPP_WITHOUT_ERRORS = [
(
"File has a single instance.",
SLOPP{Int,Int,Int}[SLOPP_INST_1],
"""
$(SLOPP_HEADER)
1
$(SLOPP_BODY_1)
"""
),
(# Two instances.
"File has two instances.",
SLOPP{Int,Int,Int}[SLOPP_INST_1, SLOPP_INST_2],
"""
$(SLOPP_HEADER)
2
$(SLOPP_BODY_1)
$(SLOPP_BODY_2)
"""
)
]

const SSSCSP_HEADER = """
***2D Rectangular Problem***
***Instances for the Single Stock Size Cutting Stock Problem (SSSCSP)***
Input parameter file: SSSCSP_parameters.txt
****************************************************************************************************
Total number of instances 
LargeObject.Length      LargeObject.Width
Number of different item types (i)
Item[i].Length  Item[i].Width   Item[i].Demand
*****************************************************************************************************"""

const SLOPP_WITH_WARNINGS = [
(
"Warning because of wrong header.",
(:warn, r"Could not match the header pattern for format:"),
"""
$(SSSCSP_HEADER)
1
$(SLOPP_BODY_1)
"""
),
(
"Warning because of extra lines after parsing all instances.",
(:warn, """The file/string indicated 1 instance(s), after consuming those some lines remained. Ignoring any further lines. The last parsed line was \"$(last(split(SLOPP_BODY_1, isequal('\n'); keepempty = false)))\" (line $(count("\n", SLOPP_HEADER) + 3 + count("\n", SLOPP_BODY_1)) considering only non-empty lines)."""),
"""
$(SLOPP_HEADER)
1
$(SLOPP_BODY_1)
Additional non-empty line that triggers warning.
"""
)
]

const SLOPP_WITH_ERRORS = [
(
"Incomplete header (not two lines of asterisks).",
GenericParseError(CPG_SLOPP{Int64,Int64,Int64}(), "Parser expects to find two lines only with asterisks; found 1 line(s) like this. As the data should be after the second of these lines, the parser did not start gathering data."),
"""
$(join(split(SLOPP_HEADER, isequal('\n'); keepempty = false)[1:end-1], "\n"))
2
$SLOPP_BODY_1
$SLOPP_BODY_2
"""
),
(
"Only the header is present (no data).",
GenericParseError(CPG_SLOPP{Int64,Int64,Int64}(), "Expected non-empty line number 10 to have 'the number of problem instances' but the file/string finished before this line."),
"""
$(SLOPP_HEADER)
"""
),
(
"Number of instances line was skipped (has wrong number of tokens).",
GenericParseError(CPG_SLOPP{Int64,Int64,Int64}(), "Expected 1 token(s) (the number of problem instances) from \"50      60\" (non-empty line number 10) got 2 token(s) instead."),
"""
$(SLOPP_HEADER)
$(SLOPP_BODY_1)
"""
),
(
"Number of instances line has an invalid token.",
GenericParseError(CPG_SLOPP{Int64,Int64,Int64}(), "Expected 'Int64' (the number of problem instances) from \"some_garbage\" (non-empty line number 10) but could not parse some of the tokens to their respective types."),
"""
$(SLOPP_HEADER)
some_garbage
$(SLOPP_BODY_1)
"""
),
(
"Provided number of instances is larger than the real number.",
GenericParseError(CPG_SLOPP{Int64,Int64,Int64}(), "Expected non-empty line number 41 to have 'L and W' but the file/string finished before this line."),
"""
$(SLOPP_HEADER)
2
$(SLOPP_BODY_1)
"""
),
(
"The original plate dimensions are not two tokens.",
GenericParseError(CPG_SLOPP{Int64,Int64,Int64}(), "Expected 2 token(s) (L and W) from \"28\" (non-empty line number 11) got 1 token(s) instead."),
"""
$(SLOPP_HEADER)
1
$(join(split(SLOPP_BODY_1, isequal('\n'); keepempty = false)[2:end], "\n"))
"""
),
(
"The original plate dimensions are two invalid tokens.",
GenericParseError(CPG_SLOPP{Int64,Int64,Int64}(), "Expected 'Int64 Int64' (L and W) from \"some garbage\" (non-empty line number 11) but could not parse some of the tokens to their respective types."),
"""
$(SLOPP_HEADER)
1
some garbage
$(join(split(SLOPP_BODY_1, isequal('\n'); keepempty = false)[2:end], "\n"))
"""
),
(
"The file ends before the number of small objects inside an instance.",
GenericParseError(CPG_SLOPP{Int64,Int64,Int64}(), "Expected non-empty line number 12 to have 'the number of small objects' but the file/string finished before this line."),
"""
$(SLOPP_HEADER)
1
100 100
"""
),
(
"Number of pieces is skipped (i.e., next line is a piece).",
GenericParseError(CPG_SLOPP{Int64,Int64,Int64}(), "Expected 1 token(s) (the number of small objects) from \"100 100\" (non-empty line number 12) got 2 token(s) instead."),
"""
$(SLOPP_HEADER)
1
100 100
100 100
"""
),
#=
( # Third line not a number.
GenericParseError(Val(:CPG_SLOPP), "dummy"),
"""
$(SLOPP_HEADER)
"""
),
( # Not enough lines for the small objects.
GenericParseError(Val(:CPG_SLOPP), "dummy"),
"""
$(SLOPP_HEADER)
"""
),
( # Not five tokens in a small object line.
GenericParseError(Val(:CPG_SLOPP), "dummy"),
"""
$(SLOPP_HEADER)
"""
),
=#
]

#= used for debug
function shallow_obj_diff(x :: T, y :: T) where {T}
	println("My test 2.")
	for field in fieldnames(T)
		@show field
		xf, yf = getfield(x, field), getfield(y, field)
		@show typeof.((xf, yf))
		if xf != yf
			for (vxf, vyf) in zip(xf, yf)
				@test vxf == vyf
			end
			println("The values in the fields $field are distinct:")
			#show(IOContext(stdout, :limit => false), xf)
			dump(IOContext(stdout, :limit => false), xf)
			println()
			#show(IOContext(stdout, :limit => false), yf)
			dump(IOContext(stdout, :limit => false), yf)
			println()
		end
	end

	return nothing
end
=#

function test_SLOPP()
	@testset "Valid SLOPP dataset, no warnings." begin
		for (case, expected, s) in SLOPP_WITHOUT_ERRORS
			# If the instance is perfectly normal then no message should be printed.
			@testset "$case" begin
				obtained = @test_logs read_from_string(Val(:CPG_SLOPP), s)
				#show(IOContext(stdout, :limit => false), obtained)
				@test expected == obtained
			end
		end
	end
	@testset "SLOPP datasets that trigger warnings." begin
		for (case, expected, s) in SLOPP_WITH_WARNINGS
			@testset "$case" begin
				@test_logs expected read_from_string(Val(:CPG_SLOPP), s)
			end
		end
	end
	@testset "SLOPP datasets that trigger exceptions." begin
		for (case, exception, s) in SLOPP_WITH_ERRORS
			@testset "$case" begin
				@test_throws exception read_from_string(Val(:CPG_SLOPP), s)
			end
		end
	end
	# TODO: tests with different number types in each field
end

const ODPW_HEADER = """
***2D Rectangular Problem***
***Problem tests for the Open Dimension Problem (ODP/W)***
Input parameter file: ODPW_parameters.txt
***********************************************************************
Total number of instances
Number of different large objects (j)
LargeObject[j].Width    LargeObject[j].Available        LargeObject[j].Value
Number of different item types (i)
Item[i].Length  Item[i].Width   Item[i].Demand
***********************************************************************"""

const MHLOPPW_HEADER = """
***2D Rectangular Problem***
***Instances for the Multiple Heterogeneous Large Object Placement Problem (MHLOPP/W)***
Input parameter file: MHLOPPW_parameters.txt
***************************************************************************************************************
Total number of instances 
Number of different large objects (j) 
LargeObject[j].Length   LargeObject[j].Width    LargeObject[j].Available
Number of different item types (i) 
Item[i].Length  Item[i].Width   Item[i].LowerBoundDemand        Item[i].UpperBoundDemand        Item[i].Value
***************************************************************************************************************"""

const SSSCSP_BODY_1 = """
82      138
20
81      7       5
6       128     2
7       127     5
39      58      1
9       130     2
28      126     3
77      87      4
6       80      1
14      80      2
35      109     2
49      62      1
5       96      1
47      6       2
32      7       2
40      109     5
33      80      5
55      114     2
10      61      5
34      138     1
61      99      5"""

const SSSCSP_BODY_2 = """
144     94
23
38      88      3
58      27      1
88      32      5
122     19      5
7       16      2
132     44      1
104     44      2
15      52      1
60      86      2
33      91      1
62      13      1
122     11      4
88      86      2
37      27      5
140     13      2
60      35      3
38      37      1
58      8       1
20      63      1
109     19      5
6       21      2
139     7       1
6       47      1"""

const SSSCSP_INST_1 = SSSCSP{Int64,Int64,Int64}(
	82, 138,
	[81, 6, 7, 39, 9, 28, 77, 6, 14, 35, 49, 5, 47, 32, 40, 33, 55,
		10, 34, 61],
	[7, 128, 127, 58, 130, 126, 87, 80, 80, 109, 62, 96, 6, 7, 109,
		80, 114, 61, 138, 99],
	[5, 2, 5, 1, 2, 3, 4, 1, 2, 2, 1, 1, 2, 2, 5, 5, 2, 5, 1, 5]
)

const SSSCSP_INST_2 = SSSCSP{Int64,Int64,Int64}(
	144, 94,
	[38, 58, 88, 122, 7, 132, 104, 15, 60, 33, 62, 122, 88, 37, 140,
		60, 38, 58, 20, 109, 6, 139, 6],
	[88, 27, 32, 19, 16, 44, 44, 52, 86, 91, 13, 11, 86, 27, 13, 35,
		37, 8, 63, 19, 21, 7, 47],
	[3, 1, 5, 5, 2, 1, 2, 1, 2, 1, 1, 4, 2, 5, 2, 3, 1, 1, 1, 5, 2, 1, 1]
)

const ODPW_BODY_1 = """
2
54      1       33
93      1       48
26
116     86      2
82      70      1
6       45      4
5       36      1
59      56      4
149     71      3
8       54      2
121     81      5
28      13      5
40      52      5
150     35      1
150     59      5
111     16      4
146     29      3
132     42      1
30      35      1
74      79      5
6       80      3
81      32      5
76      29      1
150     16      3
52      30      5
84      41      5
5       62      5
140     20      1
118     74      4"""

const ODPW_BODY_2 = """
9
133     1       104
147     1       133
71      1       36
102     1       78
145     1       108
144     1       141
87      1       76
95      1       74
57      1       54
24
142     22      1
123     20      1
30      7       2
5       58      5
24      24      2
146     53      1
117     36      2
77      21      4
142     49      5
90      22      5
57      21      4
26      42      3
150     20      1
38      52      1
111     51      5
26      42      2
19      11      5
61      21      5
57      20      5
8       22      1
48      28      1
134     40      2
67      67      5
111     26      4"""

const ODPW_INST_1 = ODPW{Int64,Int64,Int64}(
	[54, 93], [1, 1], [33, 48],
	[116, 82, 6, 5, 59, 149, 8, 121, 28, 40, 150, 150, 111, 146, 132, 30,
		74, 6, 81, 76, 150, 52, 84, 5, 140, 118],
	[86, 70, 45, 36, 56, 71, 54, 81, 13, 52, 35, 59, 16, 29, 42, 35, 79,
		80, 32, 29, 16, 30, 41, 62, 20, 74],
	[2, 1, 4, 1, 4, 3, 2, 5, 5, 5, 1, 5, 4, 3, 1, 1, 5, 3, 5, 1, 3, 5, 5,
		5, 1, 4]
)

const ODPW_INST_2 = ODPW{Int64,Int64,Int64}(
	[133, 147, 71, 102, 145, 144, 87, 95, 57],
	[1, 1, 1, 1, 1, 1, 1, 1, 1], [104, 133, 36, 78, 108, 141, 76, 74, 54],
	[142, 123, 30, 5, 24, 146, 117, 77, 142, 90, 57, 26, 150, 38, 111, 26,
		19, 61, 57, 8, 48, 134, 67, 111],
	[22, 20, 7, 58, 24, 53, 36, 21, 49, 22, 21, 42, 20, 52, 51, 42, 11,
		21, 20, 22, 28, 40, 67, 26],
	[1, 1, 2, 5, 2, 1, 2, 4, 5, 5, 4, 3, 1, 1, 5, 2, 5, 5, 5, 1, 1, 2, 5, 4]
)

const MHLOPPW_BODY_1 = """
2
138     96      39  
133     125     21  
5
5       83      0       3       415 
58      8       0       5       464 
127     9       0       3       1143
47      92      0       5       4324
77      85      0       5       6545"""

const MHLOPPW_BODY_2 = """
10
104     147     21
124     126     32
96      124     25
105     124     28
114     130     28
135     74      26
140     118     29
91      143     27
147     137     33
130     133     24
1
96      5       0       1       480"""

MHLOPPW_INST_1 = MHLOPPW{Int64,Int64,Int64}(
	[138, 133], [96, 125], [39, 21], [5, 58, 127, 47, 77],
	[83, 8, 9, 92, 85], [0, 0, 0, 0, 0], [3, 5, 3, 5, 5],
	[415, 464, 1143, 4324, 6545]
)

MHLOPPW_INST_2 = MHLOPPW{Int64,Int64,Int64}(
	[104, 124, 96, 105, 114, 135, 140, 91, 147, 130],
	[147, 126, 124, 124, 130, 74, 118, 143, 137, 133],
	[21, 32, 25, 28, 28, 26, 29, 27, 33, 24],
	[96], [5], [0], [1], [480]
)

const NON_SLOPP_NO_ERROR = Dict{Symbol, Vector{Tuple{String, String, Any}}}(
:ODPW => Tuple{String, String, Any}[
(
"File has a single instance.",
"""
$ODPW_HEADER
1
$ODPW_BODY_1
""",
ODPW{Int64,Int64,Int64}[ODPW_INST_1]
),
(
"File has two instances.",
"""
$ODPW_HEADER
2
$ODPW_BODY_1
$ODPW_BODY_2
""",
ODPW{Int64,Int64,Int64}[ODPW_INST_1, ODPW_INST_2]
),
], # END ODPW
:MHLOPPW => Tuple{String, String, Any}[
(
"File has a single instance.",
"""
$MHLOPPW_HEADER
1
$MHLOPPW_BODY_1
""",
MHLOPPW{Int64,Int64,Int64}[MHLOPPW_INST_1]
),
(
"File has two instances.",
"""
$MHLOPPW_HEADER
2
$MHLOPPW_BODY_1
$MHLOPPW_BODY_2
""",
MHLOPPW{Int64,Int64,Int64}[MHLOPPW_INST_1, MHLOPPW_INST_2]
),
], # END MHLOPPW
:SSSCSP => Tuple{String, String, Any}[
(
"File has a single instance.",
"""
$SSSCSP_HEADER
1
$SSSCSP_BODY_1
""",
SSSCSP{Int64,Int64,Int64}[SSSCSP_INST_1]
),
(
"File has two instances.",
"""
$SSSCSP_HEADER
2
$SSSCSP_BODY_1
$SSSCSP_BODY_2
""",
SSSCSP{Int64,Int64,Int64}[SSSCSP_INST_1, SSSCSP_INST_2]
),
] # END SSSCSP
) # END DICT

# Checking every possible exception for every possible type is very
# time consuming. The tests over SLOPP cover most exceptional cases.
# Here we test if the reading of valid instances of other types is
# working.
function test_non_SLOPP()
	for problem_type in (:ODPW, :MHLOPPW, :SSSCSP)
		@testset "Valid $problem_type dataset, no warnings." begin
			for (case, s, expected) in NON_SLOPP_NO_ERROR[problem_type]
				# If the instance is perfectly normal then no message should be printed.
				@testset "$case" begin
					format_symbol = Symbol("CPG_$problem_type")
					obtained = @test_logs read_from_string(Val(format_symbol), s)
					#show(IOContext(stdout, :limit => false), obtained)
					@test expected == obtained
				end
			end
		end
	end
end

test_SLOPP()
test_non_SLOPP()

