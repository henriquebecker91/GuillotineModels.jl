"""
`Data` groups the instance-parsing functions and the types involved.

The parsing is extensible, as users can define any new types and extend
the parsing functions to the types. It is advised to always return the
same type for the same format type given as input, to avoid type-unstability
and have a clear mapping.
However, "auto" formats can also be created, to simplify the user's life
(in special during testing) paying for a little overhead (and risk
of misdetection).

The parsing function is extended for parametric types (e.g.,
`CPG_SLOPP{D, S, P}`) so the data type returned (e.g., `SLOPP{D, S, P}`)
has the demand/N/M (`D`), length/width (`S`), and profit/area (`P`) fields with
the right types. The function is also extend to `Val` types for convenience.
They often map to a parametric type with sensible safe defaults (e.g.
`Val{:CPG_SLOPP}` -> `CPG_SLOPP{Int, Int, Int}`). But a finer-grained control
is also possible, as you can define something like (`Val{:CPG_SLOPP_16}` ->
`CPG_SLOPP{Int16, Int16, Int16}`). `Val` parameters are also important because
module `CommandLine` currently takes `String`s that are passed to `(Val ∘
Symbol)` in order to define the instance format to be used.

This module already extends the instance reading for many instance formats.
From the 2DCPackGen (further abbreviated here as `CPG`, Cutting and Packing
Generator, because `struct`s cannot start with a digit), the module supports:
`SLOPP`, `MHLOPPW`, `ODPW`, and `SSSCSP`. They all follow the same pattern:
`SLOPP` (the `struct` type with the data), `CPG_SLOPP{D, S, P}` (for the
clearly typed format type), and `Val{:CPG_SLOPP}` for the safe defaults (i.e.,
the same as `CPG{Int, Int, Int}`, note this will be `32` in some machined and
`64` in others). For details on these formats look at the 2DCPackGen tool.

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

export read_from_string, read_from_file

include("classic_g2kp.jl")

# TODO: restore the configuration Dict and allow:
# 	to parse instances with an id column; out of order plate ids; save ids
# TODO: SOME INSTANCES ASSUME THE READER IS COMPLETELY INSENSTIVE TO
# 	WHITESPACE, CHANGE THE METHOD TO GET THE NEXT TOKEN AND NOT WORK
# 	ORIENTED BY LINES BUT BY WORDS/TOKENS/NUMBERS
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

end # module

