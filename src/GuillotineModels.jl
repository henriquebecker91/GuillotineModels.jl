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

export build_model, get_cut_pattern, CutPattern
using DocStringExtensions # for TYPEDFIELDS
using JuMP # for JuMP.Model parameter reference

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

"""
A guillotine stage-unlimited unrestricted cutting pattern.

$(TYPEDFIELDS)

# Notes

* `cuts_are_vertical` is only relevant if `subpatterns` has length two or more.
* If `subpatterns` is non-empty, and `cuts_are_vertical == true`, then the sum
  of the `width` of the patterns in `subpatterns` cannot be greater than
  `width`, and each individual `length` of a pattern in `subpatterns` cannot be
  greater than `length`. If `cuts_are_vertical == false` the same applies
  except `length` and `width` switch roles.
* It is assumed that zero is not a valid piece index, as it is used in
  `piece_idx` to mark that a pattern is not sold as a piece (in this case,
  it is either a cut pattern or waste).
* If `piece_idx` has a value different from zero, then `subpatterns` must be
  empty. In other words, a piece cannot be subdivided.
* Waste may be represented explicitly, by having plates with `piece_idx == 0`
  but no elements in `subpatterns`, or implicitly, by having elements in
  `subpatterns` that do occupy all the area defined by `length` and `width`.
"""
struct CutPattern{D, S}
	"The length of the pattern."
	length :: S
	"The width of the pattern."
	width :: S
	"If the pattern is sold as a piece, then the piece index; otherwise zero."
	piece_idx :: D
	"The common orientation of the cuts between the plates in `subpatterns`."
	cuts_are_vertical :: Bool
	"The subpatterns that constitute the pattern."
	subpatterns :: Vector{CutPattern{D, S}}

	"Inner constructor."
	function CutPattern(
		L :: S, W :: S, piece_idx :: D, cuts_are_vertical :: Bool,
		subpatterns :: Vector{CutPattern{D, S}}
	) :: CutPattern{D, S} where {D, S}
		!iszero(piece_idx) && !isempty(subpatterns) && @error(
			"A piece cannot be subdivided; " *
			"failed: !iszero(piece_idx) && isempty(subpatterns)."
		)
		length <= zero(length) && @error(
			"The length of a pattern must be positive."
		)
		width <= zero(width) && @error(
			"The width of a pattern must be positive."
		)
		if cuts_are_vertical
			sum_width = sum(getfield.(subpatterns, :width))
			sum_width > L && @error(
				"The sum of the subpatterns width is $sum_width but the width" *
				" of the pattern is $W."
			)
			max_length = max(getfield.(subpatterns, :length))
			max_length > W && @error(
				"The largest subpattern length is $max_length but the length" *
				" of the pattern is $L."
			)
		else
			sum_length = sum(getfield.(subpatterns, :length))
			sum_length > L && @error(
				"The sum of the subpatterns length is $sum_length but the length" *
				" of the pattern is $L."
			)
			max_width = max(getfield.(subpatterns, :width))
			max_width > W && @error(
				"The largest subpattern width is $max_width but the width" *
				" of the pattern is $W."
			)
		end

		new{D, S}(L, W, piece_idx, cuts_are_vertical, subpatterns)
	end
end

"""
		CutPattern(length, width, piece_idx)

Simplified constructor for pieces and waste.
"""
function CutPattern(
	L :: S, W :: S, piece_idx :: D
) :: CutPattern{D, S} where {D, S}
	CutPattern(L, W, piece_idx, false, [])
end

"""
		CutPattern(length, width, cuts_are_vertical, subpatterns)

Simplified constructor for intermediary plates.
"""
function CutPattern(
	L :: S, W :: S, cuts_are_vertical :: Bool,
	subpatterns :: Vector{CutPattern{D, S}}
) :: CutPattern{D, S} where {D, S}
	CutPattern(L, W, zero(D), cuts_are_vertical, subpatterns)
end

"""
    get_cut_pattern(model_type, model, ::D, ::S, build_model_return)

Given a solved `model` built with `build_model(model_type, ...)` and its
respective `build_model_return`, returns a CutPattern representing the optimal
solution found. The implementation of this method is responsability of
whoever implemented the corresponding `build_model` method.

The `::Type{D}` and `::Type{S}` parameters define the integer type used for
denoting demand/'piece indexes' and the size of the pieces/plates dimensions,
the method is free to fail if those are different (especially if they are
smaller) from the `D` and `S` used for building the model.

If `build_model` returns `nothing` then `model` somehow needs to contain
all information needed to build the CutPattern. Note that the problem instance
(`l`, `w`, `L`, `W`, ...) is not passed to this method, consequently, a
common pattern is to have `build_model` return a struct containing
the problem instance and the auxiliary tables used to build the model
(especially the ones associating variables with dimensions or piece indexes).

See also: [`build_model`](@ref)
"""
function get_cut_pattern(
	model_type :: Val{T},  model :: JuMP.Model, ::Type{D}, ::Type{S},
	build_model_return :: Any
) :: CutPattern{D, S} where {T, D, S}
	@error(
		"A specialized method for $(T) should exist, but instead this " *
		" generic error fallback was called."
	)
end

include("utilities/Utilities.jl")
include("InstanceReader.jl")
#include("KnapsackPlates.jl") # not working, and no plans to fix it for now
include("flow/Flow.jl")
include("ppg2kp/PPG2KP.jl")
include("CommandLine.jl")

end # module

