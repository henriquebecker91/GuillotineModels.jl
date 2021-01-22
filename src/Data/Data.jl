"""
`Data` groups the instance-parsing functions and the types involved.

The parsing is extensible, as users can define any new types and extend
the parsing functions to the types. It is advised to always return the
same type for the same format type given as input, to avoid type-unstability
and have a clear mapping.
However, "auto" formats can also be created, to simplify the user's life
(in special during testing) paying for a little overhead (and risk
of misdetection).

The parsing function is extended for parametric types (e.g., `Classic_G2KP{D,
S, P}`) so the data type returned (e.g., `G2KP{D, S, P}`) has the demand/N/M
(`D`), length/width (`S`), and profit/area (`P`) fields with the right types.
The function is also extend to `Val` types for convenience.  They often map to
a parametric type with sensible safe defaults (e.g.  `Val{:Classic_G2KP}` ->
`Classic_G2KP{Int, Int, Int}`). But a finer-grained control is also possible,
as you can define something like (`Val{:Classic_G2KP_16}` ->
`Classic_G2KP_16{Int16, Int16, Int16}`). `Val` parameters are also important
because module `CommandLine` currently takes `String`s that are passed to `(Val
∘ Symbol)` in order to define the instance format to be used.

This module already extends the instance reading for many instance formats.
From the 2DCPackGen (further abbreviated here as `CPG`, Cutting and Packing
Generator, because `struct`s cannot start with a digit), the module supports:
`SLOPP`, `MHLOPPW`, `ODPW`, and `SSSCSP`. They all follow the same pattern:
`SLOPP` (the `struct` type with the data), `CPG_SLOPP{D, S, P}` (for the
clearly typed format type), and `Val{:CPG_SLOPP}` for the safe defaults (i.e.,
the same as `CPG{Int, Int, Int}`, note this will be `32` in some machined and
`64` in others). For details on these formats look at the 2DCPackGen tool.

Some file formats aggregate multiple instances in a single file. To allow
generic functions to deal with this possibility, it is recommended to
specialize [`GuillotineModels.Data.is_collection`](@ref). This function should
be specialized for the format type (not the return of the parsing).
A format that does not specialize it is assumed to return a single instance
upon parsing (i.e., [`GuillotineModels.Data.read_from_string`](@ref)).

The `Classic_G2KP`/`G2KP` format/data is also provided in the same fashion.
The format is described in the method documentation.

The function may `@warn` the users of "strange" (but not incorrect)
structure of the given instance, and `error` if the format is incorrect.

The current methods dispatch on the type of the first object, even if all
current format types have empty objects. Consequently, instantiations like
`CPG_SLOPP{Int8, Int16, Int32}()` and `Val(:CPG_SLOPP)` are expected,
instead of just the type itself. In the future may be possible to
add parsing options/configuration inside the format objects.
"""
module Data

import TimerOutputs.@timeit
import ..TIMER
import ..solution_value
using AutoHashEquals: @auto_hash_equals

export read_from_string, read_from_file, is_collection
export GenericParseError

# From classic_g2kp
export Classic_G2KP, G2KP
include("classic_g2kp.jl")

# From cpg
export CPG_SLOPP, CPG_MHLOPPW, CPG_ODPW, CPG_SSSCSP
export SLOPP, MHLOPPW, ODPW, SSSCSP
include("cpg.jl")

"""
    read_from_file(format, filepath :: AbstractString)

Reads the given file and apply `read_from_string` to its contents.

The (already provided) generic body of this function should be enough for most
formats without the need of explicit specialization.

See also: [`read_from_string`](@ref)
"""
@timeit TIMER function read_from_file(
	format, filepath :: AbstractString
)
	return open(filepath) do f
		read_from_string(format, read(f, String))
	end
end

"""
    is_collection(format) :: Bool

Indicate if a format is parsed to a single instance or to a collection.

If the format may store a variable number of instances, then passing the
format object to this function should return `true`; if the format always
store a single instance, then it should return false.

If this function is not specialized it assumes the format parses to a single
instance.

If the function returns `true` for some format, it is assumed that the
value returned by [`read_from_string`](@ref) (and its file counterpart)
respect the `Base.iterate` interface; if a file in this format has just a
single instance, this single instance should be wrapped in the same container
used when there are multiple instances (any overhead caused by wrapping will
probably not be larger than the overhead caused by type unstability).

Ideally, the value returned should not depend on the value of the
argument, but only the type of the argument.
"""
is_collection(format) = false

"""
    abstract type ParseError <: Exception end

Supertype for parse errors.
"""
abstract type ParseError <: Exception end

"""
    GenericParseError(format, error_message)

Most general error thrown by `parse_from_string` when it fails.

Ideally, methods extending `parse_from_string` should use a subtype of
`ParseError`, either custom or this one.
"""
struct GenericParseError <: ParseError
	format :: Any
	error_message :: String
end

function Base.showerror(io :: IO, e :: GenericParseError)
	print(io,
		"Trying to parse $(e.format) has resulted in the following error: " *
		e.error_message
	)
end

end # module

