# = 2DCPackGen =

# == Structures ==

# === Data/Problem Structures ===

@auto_hash_equals struct SLOPP{D, S, P}
	"Length of the large object."
	L :: S
	"Width of the large object."
	W :: S
	"Lengths of the small objects."
	l :: Vector{S}
	"Widths of the small objects."
	w :: Vector{S}
	"Lower bound on piece demand (minimum amount needed for a valid solution)."
	dlb :: Vector{D}
	"Upper bound on piece demand (maximum amount allowed for a valid solution)."
	dub :: Vector{D}
	"Profits of the small objects; contribution to objective function if packed."
	p :: Vector{P}
end

@auto_hash_equals struct MHLOPPW{D, S, P}
	"Lengths of the large objects."
	L :: Vector{S}
	"Widths of the large objects."
	W :: Vector{S}
	"Amount of copies available for each large object (cannot use more)."
	available :: Vector{D}
	"Lengths of the small objects."
	l :: Vector{S}
	"Widths of the small objects."
	w :: Vector{S}
	"Lower bound on piece demand (minimum amount needed for a valid solution)."
	dlb :: Vector{D}
	"Upper bound on piece demand (maximum amount allowed for a valid solution)."
	dub :: Vector{D}
	"Profit of the pieces (contribution to objective function if packed)."
	p :: Vector{P}
end

@auto_hash_equals struct ODPW{D, S, P}
	"Widths of the large objects."
	W :: Vector{S}
	"Amount of copies available for each large object (cannot use more)."
	available :: Vector{S}
	"Cost associated with the use of the large object, to be interpreted."
	cost :: Vector{S}
	"Lengths of the small objects."
	l :: Vector{S}
	"Widths of the small objects."
	w :: Vector{S}
	"The pieces to be packed (minimum amount needed for a valid solution)."
	d :: Vector{D}
end

# TODO: should this struct have the 'P' parameter?
@auto_hash_equals struct SSSCSP{D, S, P}
	"Length of the large object."
	L :: S
	"Width of the large object."
	W :: S
	"Lengths of the small objects."
	l :: Vector{S}
	"Widths of the small objects."
	w :: Vector{S}
	"The pieces to be packed (minimum amount needed for a valid solution)."
	d :: Vector{D}
end

# === Format/Instance Structures + `is_collection` ===

abstract type CPG_Format{D, S, P} end
struct CPG_SLOPP{D, S, P} <: CPG_Format{D, S, P} end
struct CPG_MHLOPPW{D, S, P} <: CPG_Format{D, S, P} end
struct CPG_ODPW{D, S, P} <: CPG_Format{D, S, P} end
struct CPG_SSSCSP{D, S, P} <: CPG_Format{D, S, P} end
is_collection(_ :: T) where {T <: CPG_Format} = true

# The Simple_CPG_SLOPP is not a real CPG format and does not return
# collections, is a simplified format that only has one instance inside
# and no header.
# TODO: in the future expand the simplified format to all supported
# CPG formats and create an abstract supertype for Simple_CPG formats.
struct Simple_CPG_SLOPP{D, S, P} end
is_collection(_ :: Val{:Simple_CPG_SLOPP}) = false
is_collection(_ :: T) where {T <: Simple_CPG_SLOPP} = false

const CPG_VALS = Union{
	Val{:CPG_SLOPP}, Val{:CPG_MHLOPPW}, Val{:CPG_ODPW}, Val{:CPG_SSSCSP}
}
is_collection(_ :: T) where {T <: CPG_VALS} = true

# === `parse_from_string` ===

# Replaces '(', ')', '*', and '/' by the escaped version.
_esc_for_regex(s) = replace(s, r"(\*|\(|\)|\/)" => s"\\\1")
_cpg_second_line(_::CPG_SLOPP) :: String = _esc_for_regex("***Instances for the Single Large Object Placement Problem (SLOPP)***")
_cpg_second_line(_::CPG_MHLOPPW) :: String = _esc_for_regex("***Instances for the Multiple Heterogeneous Large Object Placement Problem (MHLOPP/W)***")
_cpg_second_line(_::CPG_ODPW) :: String = _esc_for_regex("***Problem tests for the Open Dimension Problem (ODP/W)***")
_cpg_second_line(_::CPG_SSSCSP) :: String = _esc_for_regex("***Instances for the Single Stock Size Cutting Stock Problem (SSSCSP)***")

const CPG_DELIM = "^\\*+\$"
const CPG_DELIM_REGEX = Regex(CPG_DELIM)

function _expected_CPG_header(f) :: Regex
Regex("""
\\*\\*\\*2D Rectangular Problem\\*\\*\\*
$(_cpg_second_line(f))
Input parameter file:[^\\n]*
$CPG_DELIM
.*
$CPG_DELIM""", "mis")
end

function _check_CPG_header(f, s :: AbstractString)
	expected_pattern = _expected_CPG_header(f)
	if match(expected_pattern, s) === nothing
		@warn "Could not match the header pattern for format:" *
			"\n$(f)\n; this is:\n" * string(expected_pattern)
	end

	return nothing
end

function _find_header_end(format :: T, lines) where {T}
	qt_delims_found :: Int = 0
	current_line :: Int = 1
	last_line :: Int = length(lines)

	while qt_delims_found < 2 && current_line <= last_line
		# strip is needed because \r\n newlines keep the \r after eachline
		if match(CPG_DELIM_REGEX, strip(lines[current_line])) !== nothing
			qt_delims_found += 1
		end
		current_line += 1
	end

	if qt_delims_found < 2
		throw(GenericParseError(format,
			"Parser expects to find two lines only with asterisks;" *
			" found $(qt_delims_found) line(s) like this." *
			" As the data should be after the second of these lines," *
			" the parser did not start gathering data."
		))
	end

	return current_line
end

_cpg_eltype(::Type{CPG_SLOPP{D,S,P}}) where {D,S,P} = SLOPP{D,S,P}
_cpg_eltype(::Type{CPG_MHLOPPW{D,S,P}}) where {D,S,P} = MHLOPPW{D,S,P}
_cpg_eltype(::Type{CPG_ODPW{D,S,P}}) where {D,S,P} = ODPW{D,S,P}
_cpg_eltype(::Type{CPG_SSSCSP{D,S,P}}) where {D,S,P} = SSSCSP{D,S,P}

function _get_numbers_from_line(
	format, lines, line_number, purpose, types; check_size = false
)
	if check_size
		last_line = length(lines)
		line_number > last_line && throw(GenericParseError(format,
			"Expected non-empty line number $line_number to have '$purpose'" *
			" but the file/string finished before this line."
		))
	end

	line = lines[line_number]
	qt_numbers = length(types)
	tokens = split(line)
	qt_tokens = length(tokens)

	qt_tokens != qt_numbers && throw(GenericParseError(format,
		"Expected $qt_numbers token(s) ($purpose) from $(repr(line))" *
		" (non-empty line number $line_number) got $qt_tokens" *
		" token(s) instead."
	))

	numbers = tryparse.(types, tuple(tokens...))

	any(isequal(nothing), numbers) && throw(GenericParseError(format,
		"Expected '$(join(string.(types), ' '))' ($purpose) from $(repr(line))" *
		" (non-empty line number $line_number) but could not parse some" *
		" of the tokens to their respective types."
	))

	return numbers
end

# For both SLOPP and SSSCSP (both have a Single Large Object)
const CPG_SLO{D, S, P} = Union{CPG_SLOPP{D, S, P}, CPG_SSSCSP{D, S, P}}

# As both SLOPP and SSSCSP have the same structure of Large Object they
# can be the same method.
function _read_large_objs(
	format :: CPG_SLO{D, S, P}, lines, next :: Int, last :: Int
) :: Tuple{Int, Tuple{S, S}} where {D, S, P}
	lo_data = _get_numbers_from_line(
		format, lines, next, "L and W", (S, S); check_size = true
	) :: Tuple{S, S}

	return (next + 1, lo_data)
end

# Both ODPW and MHLOPPW have Multiple Large Objects, so we use this
# union type for the part of the code that shares that logic.
const CPG_MLO{D, S, P} = Union{CPG_ODPW{D, S, P}, CPG_MHLOPPW{D, S, P}}

function _read_large_objs(
	format :: CPG_MLO{D, S, P}, lines, next :: Int, last :: Int
) :: Tuple{Int, <: Any} where {D, S, P}
	N = (_get_numbers_from_line(
		format, lines, next, "the number of large objects", (Int,);
		check_size = true
	) :: Tuple{Int})[1] :: Int

	next += 1

	last - next + 1 < N && throw(GenericParseError(format,
		"Started parsing the large objects of an instance (at non-empty line" *
		" number $(next)) but there are not enough lines left for all" *
		" the large objects (needed $(N) there are $(last - next + 1))."
	))

	return _read_large_objs_attrs(format, lines, next, N)
end

function _read_large_objs_attrs(
	format :: CPG_ODPW{D, S, P}, lines, next :: Int, N :: Int
) :: Tuple{Int, <: Any} where {D, S, P}
	W, available, cost = (Vector{S}(undef, N), Vector{D}(undef, N),
		Vector{P}(undef, N))
	lo_vectors = (W, available, cost)

	for lo_idx = 1:N
		# No check-size needed because we checked before.
		lo_line_data = _get_numbers_from_line(
			format, lines, next, "W, available, and cost", (S, D, P)
		) :: Tuple{S, D, P}
		setindex!.(lo_vectors, lo_line_data, lo_idx)

		next += 1
	end

	return (next, lo_vectors)
end

function _read_large_objs_attrs(
	format :: CPG_MHLOPPW{D, S, P}, lines, next :: Int, N :: Int
) :: Tuple{Int, <: Any} where {D, S, P}
	L, W, available = (Vector{S}(undef, N), Vector{S}(undef, N),
		Vector{D}(undef, N))
	lo_vectors = (L, W, available)

	for lo_idx = 1:N
		# No check-size needed because we checked before.
		lo_line_data = _get_numbers_from_line(
			format, lines, next, "L, W, and available", (S, S, D)
		) :: Tuple{S, S, D}
		setindex!.(lo_vectors, lo_line_data, lo_idx)

		next += 1
	end

	return (next, lo_vectors)
end

# Both SLOPP and MHLOPPW have the same kind of Complex Single Object
# (i.e., length, width, demand lower bound, demand upper bound,
# value/profit).
const CPG_CSO{D, S, P} = Union{CPG_SLOPP{D, S, P}, CPG_MHLOPPW{D, S, P}}

# Used for both SLOPP and MHLOPPW, as both have the exact same attributes
# on an item.
function _read_small_objs_attrs(
	format :: CPG_CSO{D, S, P}, lines, next :: Int, N :: Int
) :: Tuple{Int, <: Any} where {D, S, P}
	l, w, dlb, dub, p = (Vector{S}(undef, N), Vector{S}(undef, N),
		Vector{D}(undef, N), Vector{D}(undef, N), Vector{P}(undef, N))
	so_vectors = (l, w, dlb, dub, p)

	for so_idx = 1:N
		so_line_data = _get_numbers_from_line(
			format, lines, next, "length, width, demand LB, demand UB, and value",
			(S, S, D, D, P) # No check-size needed because we checked before.
		) :: Tuple{S, S, D, D, P}
		setindex!.(so_vectors, so_line_data, so_idx)

		next += 1
	end

	return (next, (l, w, dlb, dub, p))
end

# Both ODPW and SSSCSP have a Simple Small Objects (just length, width,
# and demand). So they can share the extraction of these attributes.
const CPG_SSO{D, S, P} = Union{CPG_ODPW{D, S, P}, CPG_SSSCSP{D, S, P}}

function _read_small_objs_attrs(
	format :: CPG_SSO{D, S, P}, lines, next :: Int, N :: Int
) :: Tuple{Int, <: Any} where {D, S, P}
	l, w, d = (Vector{S}(undef, N), Vector{S}(undef, N),
		Vector{D}(undef, N))
	so_vectors = (l, w, d)

	for so_idx = 1:N
		so_line_data = _get_numbers_from_line(
			format, lines, next, "length, width, demand",
			(S, S, D) # No check-size needed because we checked before.
		) :: Tuple{S, S, D}
		setindex!.(so_vectors, so_line_data, so_idx)

		next += 1
	end

	return (next, (l, w, d))
end

function _read_small_objs(
	format :: CPG_Format{D, S, P}, lines, next :: Int, last :: Int
) :: Tuple{Int, <: Any} where {D, S, P}
	N = (_get_numbers_from_line(
		format, lines, next, "the number of small objects", (Int,);
		check_size = true
	) :: Tuple{Int})[1] :: Int

	next += 1

	last - next + 1 < N && throw(GenericParseError(format,
		"Started parsing the small objects of an instance (at non-empty line" *
		" number $(next)) but there are not enough lines left for all" *
		" the small objects (needed $(N) there are $(last - next + 1))."
	))

	return _read_small_objs_attrs(format, lines, next, N)
end

function _read_instance(
	format :: T, lines, next :: Int, last :: Int
) where {T <: CPG_Format}
	next, lo = _read_large_objs(format, lines, next, last)
	next, so = _read_small_objs(format, lines, next, last)
	return next, _cpg_eltype(T)(lo..., so...)
end

# The format is passed to build error objects if necessary, otherwise it
# is not used.
function _read_number_of_instances(
	@nospecialize(format :: CPG_Format), lines, next :: Int
) :: Int # returns number of instances
	qt_instances = (_get_numbers_from_line(
		format, lines, next, "the number of problem instances", (Int,);
		check_size = true
	) :: Tuple{Int})[1] :: Int

	return qt_instances
end

"""
    read_from_string(format :: CPG_Format{D, S, P}, s :: AbstractString)

Expect a string in the format given and parse it to the adequate instance
object. The `format` may be `CPG_SLOPP{D, S, P}()`, `CPG_MHLOPPW{D, S, P}()`,
`CPG_ODPW{D, S, P}()`, `CPG_SSSCSP{D, S, P}()`. All these formats come from
2DCPackGen (10.1016/j.ejor.2014.02.059, see
https://sites.google.com/gcloud.fe.up.pt/cutting-and-packing-tools/2dcpackgen
for the code). Each of these formats starts with a self-explanatory header.
These headers are copyied below.

```
***2D Rectangular Problem***
***Instances for the Single Large Object Placement Problem (SLOPP)***
Input parameter file: SLOPP_parameters.txt
***************************************************************************************************************
Total number of instances
LargeObject.Length      LargeObject.Width
Number of different item types (i)
Item[i].Length  Item[i].Width   Item[i].LowerBoundDemand        Item[i].UpperBoundDemand        Item[i].Value
***************************************************************************************************************
```

```
***2D Rectangular Problem***
***Instances for the Multiple Heterogeneous Large Object Placement Problem (MHLOPP/W)***
Input parameter file: MHLOPPW_parameters.txt
***************************************************************************************************************
Total number of instances 
Number of different large objects (j) 
LargeObject[j].Length   LargeObject[j].Width    LargeObject[j].Available
Number of different item types (i) 
Item[i].Length  Item[i].Width   Item[i].LowerBoundDemand        Item[i].UpperBoundDemand        Item[i].Value
***************************************************************************************************************
```

```
***2D Rectangular Problem***
***Problem tests for the Open Dimension Problem (ODP/W)***
Input parameter file: ODPW_parameters.txt
***********************************************************************
Total number of instances 
Number of different large objects (j)
LargeObject[j].Width    LargeObject[j].Available        LargeObject[j].Value
Number of different item types (i)
Item[i].Length  Item[i].Width   Item[i].Demand
***********************************************************************
```

```
***2D Rectangular Problem***
***Instances for the Single Stock Size Cutting Stock Problem (SSSCSP)***
Input parameter file: SSSCSP_parameters.txt
****************************************************************************************************
Total number of instances
LargeObject.Length      LargeObject.Width
Number of different item types (i)
Item[i].Length  Item[i].Width   Item[i].Demand
*****************************************************************************************************
```

The `D`, `S`, and `P` type parameters indicate the integer type which
should be used to store `D`emand, `S`ize, and `P`rofit (or piece/plate area).

Returns an object of type `YYY{D, S, P}` in which `YYY` is what comes after
`CPG_` in the given `format` type.
"""
function read_from_string(
	format :: CPG_Format{D, S, P}, s :: AbstractString
) where {D, S, P}
	_check_CPG_header(format, s)
	lines = split(s, isequal('\n'); keepempty = false)
	last_line_number = length(lines)
	line_number_after_header = _find_header_end(format, lines)

	qt_instances = _read_number_of_instances(
		format, lines, line_number_after_header
	)
	curr_line_number = line_number_after_header + 1

  instances = Vector{_cpg_eltype(typeof(format))}(undef, qt_instances)

	for curr_inst_number = 1:qt_instances
		curr_line_number, instances[curr_inst_number] =
			_read_instance(format, lines, curr_line_number, last_line_number) ::
				Tuple{Int, _cpg_eltype(typeof(format))}
	end

	if curr_line_number <= last_line_number
		@warn "The file/string indicated $(qt_instances) instance(s), after" *
			" consuming those some lines remained. Ignoring any further lines." *
			" The last parsed line was \"$(lines[curr_line_number - 1])\"" *
			" (line $(curr_line_number - 1) considering only non-empty lines)."
	end

	return instances
end

"""
    read_from_string(_ :: Simple_CPG_SLOPP{D, S, P}, s :: AbstractString)

Similar to `read_from_string(CPG_SLOPP{D, S, P}(), s)` which is described
in `read_from_string(format :: CPG_Format{D, S, P}, s :: AbstractString)`,
the only difference in format compared to `CPG_SLOPP` is that
`Simple_CPG_SLOPP` format expects just a single instance and no header.
"""
function read_from_string(
	format :: Simple_CPG_SLOPP{D, S, P}, s :: AbstractString
) where {D, S, P}
	lines = split(s, isequal('\n'); keepempty = false)
	last_line_number = length(lines)

	curr_line_number, instance =
		_read_instance(CPG_SLOPP{D, S, P}(), lines, 1, last_line_number) ::
			Tuple{Int, _cpg_eltype(CPG_SLOPP{D, S, P})}

	if curr_line_number <= last_line_number
		@warn "After consuming all lines of the instance" *
			" some lines remained. Ignoring any further lines." *
			" The last parsed line was \"$(lines[curr_line_number - 1])\"" *
			" (line $(curr_line_number - 1) considering only non-empty lines)."
	end

	return instance
end

"""
    read_from_string(::Val{:Simple_CPG_SLOPP}, s :: AbstractString)

Convenience method. To be called as `read_from_string(Val(:Simple_CPG_SLOPP), s)`.
Equivalent to the call `read_from_string(Simple_CPG_SLOPP{Int, Int, Int}(), s)`.
"""
function read_from_string(
	_ :: Val{:Simple_CPG_SLOPP}, s :: AbstractString
)
	return read_from_string(Simple_CPG_SLOPP{Int, Int, Int}(), s)
end

"""
    read_from_string(::Val{:CPG_SLOPP}, s :: AbstractString)

Convenience method. To be called as `read_from_string(Val(:CPG_SLOPP), s)`.
Equivalent to the call `read_from_string(CPG_SLOPP{Int, Int, Int}(), s)`.
"""
function read_from_string(
	_ :: Val{:CPG_SLOPP}, s :: AbstractString
)
	return read_from_string(CPG_SLOPP{Int, Int, Int}(), s)
end

"""
    read_from_string(::Val{:CPG_ODPW}, s :: AbstractString)

Convenience method. To be called as `read_from_string(Val(:CPG_ODPW), s)`.
Equivalent to the call `read_from_string(CPG_ODPW{Int, Int, Int}(), s)`.
"""
function read_from_string(
	_ :: Val{:CPG_ODPW}, s :: AbstractString
)
	return read_from_string(CPG_ODPW{Int, Int, Int}(), s)
end

"""
    read_from_string(::Val{:CPG_MHLOPPW}, s :: AbstractString)

Convenience method. To be called as `read_from_string(Val(:CPG_MHLOPPW), s)`.
Equivalent to the call `read_from_string(CPG_MHLOPPW{Int, Int, Int}(), s)`.
"""
function read_from_string(
	_ :: Val{:CPG_MHLOPPW}, s :: AbstractString
)
	return read_from_string(CPG_MHLOPPW{Int, Int, Int}(), s)
end

"""
    read_from_string(::Val{:CPG_SSSCSP}, s :: AbstractString)

Convenience method. To be called as `read_from_string(Val(:CPG_SSSCSP), s)`.
Equivalent to the call `read_from_string(CPG_SSSCSP{Int, Int, Int}(), s)`.
"""
function read_from_string(
	_ :: Val{:CPG_SSSCSP}, s :: AbstractString
)
	return read_from_string(CPG_SSSCSP{Int, Int, Int}(), s)
end

# === `write_to_file` implementations ===

function write_to_file(
	format :: Val{:Simple_CPG_SLOPP},
	instance :: SLOPP{D, S, P},
	io :: IO
) where {D, S, P}
	return write_to_file(Simple_CPG_SLOPP{Int, Int, Int}(), instance, io)
end

@timeit TIMER function write_to_file(
	format :: Simple_CPG_SLOPP{D, S, P},
	instance :: SLOPP{D, S, P},
	io :: IO
) where {D, S, P}
	@unpack L, W, l, w, dlb, dub, p = instance
	write(io, "$L $W\n")
	write(io, "$(length(p))\n")
	for (pl, pw, pdlb, pdub, pp) in zip(l, w, dlb, dub, p)
		write(io, "$pl $pw $pdlb $pdub $pp\n")
	end
	flush(io)
	return io
end

function write_to_file(
	format :: Val{:CPG_MHLOPPW}, instance, io :: IO
) where {D, S, P}
	return write_to_file(CPG_MHLOPPW{Int, Int, Int}(), instance, io)
end

const MHLOPPW_HEADER = """
***2D Rectangular Problem***
***Instances for the Multiple Heterogeneous Large Object Placement Problem (MHLOPP/W)***
Input parameter file: NONE
***************************************************************************************************************
Total number of instances 
Number of different large objects (j)
LargeObject[j].Length	LargeObject[j].Width	LargeObject[j].Available
Number of different item types (i)
Item[i].Length	Item[i].Width	Item[i].LowerBoundDemand	Item[i].UpperBoundDemand	Item[i].Value
***************************************************************************************************************"""

@timeit TIMER function write_to_file(
	format :: CPG_MHLOPPW{D, S, P},
	instance :: MHLOPPW{D, S, P},
	io :: IO
) where {D, S, P}
	# Ideally we should wrap on a StaticArray because there is no reason for this
	# allocation, but this is not performance-critical code.
	write_to_file(CPG_MHLOPPW{D, S, P}(), [instance], io)
	return io
end

@timeit TIMER function write_to_file(
	format :: CPG_MHLOPPW{D, S, P},
	instances :: AbstractVector{MHLOPPW{D, S, P}},
	io :: IO
) where {D, S, P}
	println(io, MHLOPPW_HEADER)
	write(io, "$(length(instances))\n")
	for instance in instances
		@unpack L, W, available, l, w, dlb, dub, p = instance
		qt_plates = only(unique(length.((L, W, available))))
		write(io, "$qt_plates\n")
		for (PL, PW, PD) in zip(L, W, available)
			write(io, "$PL $PW $PD\n")
		end
		qt_pieces = only(unique(length.((l, w, dlb, dub, p))))
		write(io, "$qt_pieces\n")
		for (pl, pw, pdlb, pdub, pp) in zip(l, w, dlb, dub, p)
			write(io, "$pl $pw $pdlb $pdub $pp\n")
		end
	end
	flush(io)
	return io
end

# CPG_SSSCSP write_to_file

function write_to_file(
	format :: Val{:CPG_SSSCSP}, instance, io :: IO
) where {D, S, P}
	return write_to_file(CPG_SSSCSP{Int, Int, Int}(), instance, io)
end

const SSSCSP_HEADER = """
***2D Rectangular Problem***
***Instances for the Single Stock Size Cutting Stock Problem (SSSCSP)***
Input parameter file: NONE
****************************************************************************************************
Total number of instances 
LargeObject.Length	LargeObject.Width
Number of different item types (i)
Item[i].Length	Item[i].Width	Item[i].Demand
*****************************************************************************************************"""

@timeit TIMER function write_to_file(
	format :: CPG_SSSCSP{D, S, P},
	instance :: SSSCSP{D, S, P},
	io :: IO
) where {D, S, P}
	# Ideally we should wrap on a StaticArray because there is no reason for this
	# allocation, but this is not performance-critical code.
	write_to_file(CPG_SSSCSP{D, S, P}(), [instance], io)
	return io
end

@timeit TIMER function write_to_file(
	format :: CPG_SSSCSP{D, S, P},
	instances :: AbstractVector{SSSCSP{D, S, P}},
	io :: IO
) where {D, S, P}
	write(io, SSSCSP_HEADER)
	write(io, '\n')
	write(io, "$(length(instances))\n")
	for instance in instances
		@unpack L, W, l, w, d = instance
		write(io, "$(L)\t$(W)\n")
		write(io, "$(length(d))\n")
		for (li, wi, di) in zip(l, w, d)
			write(io, "$(li)\t$(wi)\t$(di)\n")
		end
	end
	flush(io)
	return io
end

# === `GuillotineModels.solution_value` implementations

import ..CutPattern
function solution_value(
	:: Val{:G2KP}, i :: SLOPP{D, S, P}, cp :: CutPattern{D, S}
) where {D, S, P}
	@assert iszero(sum(i.dlb))
	g2kp_instance = G2KP{D, S, P}(i.L, i.W, i.l, i.w, i.dub, i.p)
	return solution_value(Val(:G2KP), g2kp_instance, cp)
end

# G2OPP does not have a solution value, because it is not an optimization
# problem but instead a decision problem. Ideally, a formulation that solves
# it will have no objective function set (or better, all coefficients of the
# objective function will bet set to the default, that is zero) and therefore
# we always return zero to match what the solver will return.
function solution_value(
	:: Val{:G2OPP}, instance, cp :: CutPattern{D, S}
) where {D, S}
	return 0
end

# The instance is taken to adhere to the signature but it is unneeded
# to compute the solution value.
function solution_value(
	:: Val{:G2CSP}, instance, bins :: Vector{CutPattern{D, S}}
) where {D, S}
	return length(bins)
end

function solution_value(
	:: Val{:G2MKP}, i :: MHLOPPW{D, S, P}, cps :: Vector{CutPattern{D, S}}
) where {D, S, P}
	@assert iszero(sum(i.dlb))
	g2kp_instance = G2KP{D, S, P}(only(i.L), only(i.W), i.l, i.w, i.dub, i.p)
	value = zero(P)
	for cp in cps
		value += solution_value(Val(:G2KP), g2kp_instance :: G2KP{D, S, P}, cp)
	end

	return value
end

