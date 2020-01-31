"""
GuillotineModels is a collection of mathematical models for 2D Guillotine
Cutting problems, all implemented using JuMP. It was developed as part of
the PhD thesis of Henrique Becker.

The main features are:
...
"""
module GuillotineModels

# Just include all submodules. Except by the specializations of build_model,
# all methods are part of some submodule and not of the top-level.
# The module names are not exported, it is expected the user will abbreviate
# the module used with a constant global variable.

export build_model

"""
"""
function build_model(
	::Val{T}, model, d :: Vector{D}, p :: Vector{P},
	l :: Vector{S}, w :: Vector{S}, L :: S, W :: S,
	options :: Dict{String, Any} = Dict{String, Any}()
) where {T, D, S, P}
	@error(
		"A specialized method for $(T) should exist, but instead this " *
		" generic error fallback was called."
	)
end

include("utilities/Utilities.jl")
include("InstanceReader.jl")
include("KnapsackPlates.jl")
include("flow/Flow.jl")
include("ppg2kp/PPG2KP.jl")
include("CommandLine.jl")

end # module

