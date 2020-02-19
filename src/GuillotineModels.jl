# TODO: for now, many documentation strings do not wrap at 79 characters
# because this break list items, discover where is the problem (incorrect use
# of markdown? documenter.jl is broken?).
"""
GuillotineModels is a collection of mathematical models for 2D Guillotine
Cutting problems, all implemented using Julia+JuMP.

It was developed as part of the PhD thesis of [Henrique Becker]
(https://orcid.org/0000-0003-3879-2691).

The main features are:
1. The implementation of distinct models using the same technology (Julia+JuMP) which is solver-agnostic.
2. The implementation of a CommandLine interface that make it easy to call the implemented models from the command-line to be solved by an specified solver, and also is extendable for new models and solvers.
"""
module GuillotineModels

# Just include all submodules. Except by the specializations of build_model,
# all methods are part of some submodule and not of the top-level.

export build_model

"""
    build_model(::Val{T}, model, d, p, l, w, L, W[, options])

Given an instance of the demand-constrained unrestricted stage-unlimited
knapsack 2D cutting problem without rotation (defined by `d`, `p`, `l`, `w`,
`L`, `W`), add variables and constraints to `model` to make it a `::Val{T}`
model for the specified instance. A dict of options that are
`::Val{T}`-specific may be also provided (otherwise defaults are used).
The values returned are `::Val{T}`-specific and are often byproducts of
the model construction that are needed to assemble a solution from the
variables values of a solved model.

# Arguments

1. `::Val{T}`: Must be a `value-type` of a symbol identifying the mathematical to be built.
2. `model`: Something that behaves as a JuMP model, it wil have variables and constraints added to it.
3. `d::Vector{D}`: The demand of the pieces.
4. `p::Vector{D}`: The profit of the pieces.
5. `l::Vector{D}`: The length of the pieces.
6. `w::Vector{D}`: The width of the pieces.
7. `L::S`: The length of the original plate.
8. `W::S`: The width of the original plate.
9. `options::Dict{String, Any}`: Model-specific options.
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

