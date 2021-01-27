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

import Dates
import TimerOutputs
const TIMER = TimerOutputs.TimerOutput()

# Just include all submodules. Except by the specializations of build_model,
# all methods are part of some submodule and not of the top-level.

export build_model, get_cut_pattern, CutPattern, to_pretty_str, simplify!
export cut_pattern_profit
export TimeoutError
using DocStringExtensions # for TYPEDFIELDS
using JuMP # for JuMP.Model parameter reference

"""
    TimeoutError <: Exception

An exception subtype representing a time limit not respected.
"""
struct TimeoutError <: Exception
	"The moment the timing started (in seconds since epoch)."
	start :: Float64
	"The amount of seconds the procedure had to finish."
	limit :: Float64
	"The moment the code noticed the timeout (in seconds since epoch)"
	now :: Float64
end

function Base.showerror(io :: IO, e :: TimeoutError)
	start, limit, now = round.((e.start, e.limit, e.now); digits = 3)
	since_start = round(now - start; digits = 3)
	since_limit = round(since_start - limit; digits = 3)
	print(io,
		"TimeoutError: A timer started at $(Dates.unix2datetime(start)) with" *
		" limit of $limit seconds (i.e., until $(Dates.unix2datetime(start +
		limit))) was discovered running at $(Dates.unix2datetime(now)), this is," *
		" $(since_start) seconds since start or $(since_limit) seconds" *
		" since it should have stopped.\n"
	)
end

"""
    throw_if_timeout(start, limit, now)

Throws a TimeoutError if `start - now > limit`.
"""
function throw_if_timeout(start, limit, now)
	(now - start) > limit && throw(TimeoutError(start, limit, now))
	return
end

"""
    throw_if_timeout_now(start, limit)

Throws a TimeoutError if `start - time() > limit`.
"""
function throw_if_timeout_now(start, limit)
	throw_if_timeout(start, limit, time())
	return
end

@enum BuildStopReason BUILT_MODEL FOUND_OPTIMUM

"""
    build_model(problem, formulation, instance, model[, options])

Given an `instance` of `problem`, build `formulation` inside `model`.

A dict of options that are `problem`- and `formulation`-specific may
be also provided (otherwise defaults are used).

Should always return two values: the first is always an element of
`BuildStopReason`; the second is `problem`- and `formulation`-specific.

# Arguments

1. `problem`: An object that identifies the problem being solved, e.g., `Val(:G2KP)`. Each problem often requires a different `instance` type but sometimes an instance of a more general problem (that is not using the extra generality) may be passed to a more specific problem.
2. `formulation`: An object that identifies the formulation used to model the `problem`, e.g., `Val(:PPG2KP)`. Often same formulation can be adapted to multiple different `problem`s with minimal changes.
3. `instance`: An object containing the `problem` data, e.g., `G2KP` or `SLOPP`.
4. `model`: An object that behaves as a JuMP model, it will have variables and constraints added to it.
5. `options::Dict{String, Any}`: Problem- and formulation-specific options.

# Returns

1. `::BuildStopReason`: The reason for which the method returned. It can be: BUILT_MODEL (the model was built successfully) or FOUND_OPTIMUM (in the process of building the model, the optimal solution was already found, so the building process was abandoned). See [`get_cut_pattern`](@ref) for more information.
2. `::Any`: the specific type of the second returned value depends on the `problem` and `formulation`, however it should always exist (even if it is `nothing`). It is passed as an argument to `get_cut_pattern` in order to assemble a solution (CutPattern) from the variables values of a model.
"""
function build_model end

#=
The list item below was removed because for now there seems to be no value
in having a field indicating the piece was rotated or not.
 * `piece_is_rotated` is only relevant if `piece_idx` is not zero; it
  should always be `false` if `piece_idx` is zero. If `true` then the piece
	dimensions are considered swapped in the context of the pattern.
=#
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
  `subpatterns` that do not fill all the area defined by `length` and `width`.
"""
struct CutPattern{D, S}
	"The length of the pattern."
	length :: S
	"The width of the pattern."
	width :: S
	"If the pattern represents a piece, then the piece index; otherwise zero."
	piece_idx :: D
	#"If the piece the pattern represents was rotated (90 degrees)."
	#piece_is_rotated :: Bool
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
		L <= zero(L) && @error(
			"The length of a pattern must be positive."
		)
		W <= zero(W) && @error(
			"The width of a pattern must be positive."
		)
		if cuts_are_vertical
			sum_width = reduce(+, getfield.(subpatterns, :width); init = zero(S))
			sum_width > W && @error(
				"The sum of the subpatterns width is $sum_width but the width" *
				" of the pattern is $W."
			)
			max_length = reduce(max, getfield.(subpatterns, :length); init = zero(S))
			max_length > L && @error(
				"The largest subpattern length is $max_length but the length" *
				" of the pattern is $L."
			)
		else
			sum_length = reduce(+, getfield.(subpatterns, :length); init = zero(S))
			sum_length > L && @error(
				"The sum of the subpatterns length is $sum_length but the length" *
				" of the pattern is $L."
			)
			max_width = reduce(max, getfield.(subpatterns, :width); init = zero(S))
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
	L :: S, W :: S, piece_idx :: D = 0
) :: CutPattern{D, S} where {D, S}
	CutPattern(L, W, piece_idx, false, CutPattern{D, S}[])
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
    solution_value(problem, instance, solution)

Return the value of `solution` for `problem` and `instance`.

The returned value is often an `Int`. If the solution (i) was extracted
from a model with a valid primal solution and (ii) the model correctly
computed the objective function, then the value returned by this function
should match the value of the objective function in the model from which
the solution was extracted. Basically, this function serves as a
formulation-agnostic double-check of the obtained solution value.

`solution` is something returned by `get_cut_pattern` so often it is
either a `Vector{CutPattern{D, S}}` or a single `CutPattern{D, S}`.

!!! This function does not check the validity of the given solution,
it assumes a valid solution and just compute the value it would have
if it was valid.
"""
function solution_value end

# INTERNAL USE. See `remove_waste!`.
function inner_rec_remove_waste!(p :: CutPattern{D, S}) :: Bool where {D, S}
	filter!(inner_rec_remove_waste!, p.subpatterns)
	return !iszero(p.piece_idx) || !isempty(p.subpatterns)
end

# INTERNAL USE.
# Removal every non-root pattern that has no piece among its descendants.
function remove_waste!(p :: CutPattern{D, S}) :: CutPattern{D, S} where {D, S}
	inner_rec_remove_waste!(p)
	return p
end

# INTERNAL USE.
# Remove every non-root pattern that has no piece among its descendants.
function promote_single_childs!(
	p :: CutPattern{D, S}
) :: CutPattern{D, S} where {D, S}
	for (index, value) in pairs(p.subpatterns)
		p.subpatterns[index] = promote_single_childs!(value)
	end
	isone(length(p.subpatterns)) && return first(p.subpatterns)
	return p
end

# INTERNAL USE.
# For every subpattern that have itself subpatterns and shares the value of
# `cuts_are_vertical` with its parent pattern, replace the subpattern with
# their immediate children.
function flatten!(p :: CutPattern{D, S}) :: CutPattern{D, S} where {D, S}
	isempty(p.subpatterns) && return p
	# Check beforehand if there is any children to be flattened save it.
	to_flatten = @. (
		(p.cuts_are_vertical === getfield(p.subpatterns, :cuts_are_vertical)) &
		iszero(getfield(p.subpatterns, :piece_idx))
	)
	need_flattening = any(to_flatten)
	need_flattening && (aux_vector = Vector{CutPattern{D, S}}())
	for (i, subpatt) in pairs(p.subpatterns)
		flatten!(subpatt)
		if need_flattening
			if to_flatten[i]
				append!(aux_vector, subpatt.subpatterns)
			else
				push!(aux_vector, subpatt)
			end
		end
	end
	if need_flattening
		# The CutPattern object is immutable so we cannot just attribute the new
		# vector to `subpatterns` field. The trick below becomes necessary.
		empty!(p.subpatterns)
		append!(p.subpatterns, aux_vector)
	end
	return p
end

"""
    simplify!(p :: CutPattern{D, S}) :: CutPattern{D, S} where {D, S}

This method is free to make any changes while keeping all the packed pieces and
the the validity of the solution, its purpose is to make the structure easier
to understand by human beings. If you are trying to debug a model you will
probably want to see both the simplified and not simplified output. The
simplified will make it easier to understand the solution itself, and the
non-simplified will show peculiarities that given insight on how the model
works/'see things' internally.

NOTE: the returned pattern is not always the given parameter, in some rare
cases (a root plate with a single piece inside) the pattern returned may
be another object.

NOTE: this method does not create any new `CutPattern` objects but instead
prunes and rearranges the pattern tree (i.e., throws away some objects and
changes the position of others in the tree).
"""
function simplify!(p :: CutPattern{D, S}) :: CutPattern{D, S} where {D, S}
	# The order these subroutines are called is relevant. For maximum
	# simplification we need to remove waste before promoting single childs
	# (or some childs that will have no siblings after waste removal will not
	# be promoted), and to promote single childs before flattening (or an
	# alternating pattern H-V-H-'2 pieces' will end up like H-H-'2 pieces'
	# and not like H-'2 pieces' as it should).
	p = remove_waste!(p)
	p = promote_single_childs!(p)
	p = flatten!(p)
	return p
end

"""
    to_pretty_str(p :: CutPattern{D, S}; kwargs...) :: String

Creates a simplified and indented string representing the CutPattern.

# Format

The format represents pieces as "`piece_idx`p`length`x`width`" (e.g., "1p10x20"
is a copy of piece 1 which have length 10 and width 20); non-piece patterns are
represented by "P`length`x`width`" (note the `P` is uppercase); if the
non-piece pattern has subpatterns (i.e., is not waste) then it starts a set of
vertical (horizontal) cuts with `[` (`{`) and close it with `]` (`}`). There is
always whitespace between the elements of such sets but, for conciseness and
ease of reading, if all the elements of a subpattern have no children they are
separated by single spaces (no matter how long the list), otherwise they are
separated by newlines.

# Keyword Arguments

* `lvl :: Int64`: The current level of indentation.
* `indent_str`: The string or character to be repeated as indentation.

"""
function to_pretty_str(
	p :: CutPattern{D, S};
	lvl = 0,
	indent_str = "  " # i.e., two spaces
) :: String where {D, S}
	indent = indent_str^lvl
	!iszero(p.piece_idx) && return "$indent$(p.piece_idx)p$(p.length)x$(p.width)"
	ob = p.cuts_are_vertical ? '[' : '{'
	cb = p.cuts_are_vertical ? ']' : '}'
	base = "$(indent)P$(p.length)x$(p.width)"
	isempty(p.subpatterns) && return base
	if any(!isempty, getfield.(p.subpatterns, :subpatterns))
		joined_sub = mapreduce(*, p.subpatterns) do subpatt
			to_pretty_str(subpatt; lvl = lvl + 1, indent_str = indent_str) * '\n'
		end
		return "$base$ob\n$joined_sub$indent$cb"
	else
		piece_strs = map(p.subpatterns) do piece
			to_pretty_str(piece; lvl = 0)
		end
		joined_pieces = join(piece_strs, " ")
		return "$base$ob $joined_pieces $cb"
	end
end

"""
    get_cut_pattern(problem, formulation, model, build_model_return)

Given a `model` built with `build_model(problem, formulation, ...)` and its
respective `build_model_return` (i.e., the second return of `build_model`),
returns a CutPattern representing the optimal solution found. If the first
return of `build_model` is `BUILT_MODEL`, then the model needs to be solved
before calling this method; however, if it is `FOUND_OPTIMUM` then the model
does not need to be solved.

The implementation of this method is responsability of
whoever implemented the corresponding `build_model` method.

The `D` and `S` type parameters from `CutPattern{D, S}` are inferred from
`build_model_return` as of now.

If `build_model` returns `nothing` then `model` somehow needs to contain
all information needed to build the CutPattern. Note that the problem instance
(`l`, `w`, `L`, `W`, ...) is not passed to this method, consequently, a
common pattern is to have `build_model` return a struct containing
the problem instance and the auxiliary tables used to build the model
(especially the ones associating variables with dimensions or piece indexes).

The solutions of some problems are not a single `CutPattern` but instead
multiple `CutPatterns` inside an iterable collection of them (often a
`Vector`).

See also: [`build_model`](@ref)
"""
function get_cut_pattern end

include("Utilities/Utilities.jl")
include("Data/Data.jl")
#include("KnapsackPlates.jl") # not working, and no plans to fix it for now
include("Flow/Flow.jl")
include("PPG2KP/PPG2KP.jl")
include("CommandLine.jl")

end # module
