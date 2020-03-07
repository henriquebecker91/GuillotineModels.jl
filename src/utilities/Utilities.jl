module Utilities

using DocStringExtensions # for TYPEDFIELDS

"""
Grouping of the length and width piece vectors (in original and sorted order),
and their reverse indexes, allowing to, for example, iterate the pieces
by length while having `O(1)` access to their width.

$(TYPEDFIELDS)

"""
struct SortedLinkedLW{D, S}
	"The pieces length in the original order."
	l :: Vector{S}
	"The pieces width in the original order."
	w :: Vector{S}
	"The pieces length sorted by increase-or-same order."
	sl :: Vector{S}
	"The pieces width sorted by increase-or-same order."
	sw :: Vector{S}
	"Translator from indexes in `sl` to piece index (`l` and `w`)."
	sli2pii :: Vector{D}
	"Translator from indexes in `sw` to piece index (`l` and `w`)."
	swi2pii :: Vector{D}
	"Translator from piece indexes (`l` and `w`) to index in `sl`."
	pii2sli :: Vector{D}
	"Translator from piece indexes (`l` and `w`) to index in `sw`."
	pii2swi :: Vector{D}
end
export SortedLinkedLW

"""
    SortedLinkedLW(::Type{D}, l :: [S], w :: [S])

Construts a SortedLinkedLW structure using type `D` as the type for indexes,
and `l` and `w` as the pieces length and width in the 'original' ordering.

NOTE: `l` and `w` are not copyed, so mutating them after  will silently and
completely invalidate the entire structure. They are not copyed for performance
reasons and because if the original vectors may be changed then the concept of
'original order' is not really relevant.
"""
function SortedLinkedLW(::Type{D}, l :: Vector{S}, w :: Vector{S}) where {D, S}
	@assert length(l) == length(w)
	n = length(l)
	sl = sort(l)
	sw = sort(w)
	sli2pii = sort!(collect(1:n), by = pii -> l[pii])
	swi2pii = sort!(collect(1:n), by = pii -> w[pii])
	pii2sli = Vector{D}(undef, n)
	pii2swi = Vector{D}(undef, n)
	for si = 1:n
		pii2sli[sli2pii[si]] = si
		pii2swi[swi2pii[si]] = si
	end
	SortedLinkedLW(l, w, sl, sw, sli2pii, swi2pii, pii2sli, pii2swi)
end

using JuMP
using TimerOutputs

# Style guideline: as the module block is left unindented, the @timeit
# blocks that wrap the whole method body also are not indented.

"""
    num_all_constraints(m) :: Int64

JuMP only allow to query the number of constraints of some specific type;
this method queries all constraint types used in the model and then sums the
number of constraints of each type.
"""
function num_all_constraints(m) :: Int64
	sum = 0 :: Int64
	for (ftype, stype) in list_of_constraint_types(m)
		sum += num_constraints(m, ftype, stype)
	end
	return sum
end
export num_all_constraints

# see https://github.com/JuliaOpt/MathOptInterface.jl/issues/776
function reduced_cost(var) :: Float64
	rc = 0.0
	has_upper_bound(var) && (rc += shadow_price(UpperBoundRef(var)))
	has_lower_bound(var) && (rc += shadow_price(LowerBoundRef(var)))
	is_fixed(var) && (rc += shadow_price(FixRef(var)))
	!has_upper_bound(var) && !has_lower_bound(var) && !is_fixed(var) &&
		@warn "CAUTION: reduce_cost was called over a variable with no bounds"
	rc
end
export reduced_cost

"""
Stores the type and bounds of a variable so they may be restored.

$(TYPEDFIELDS)

"""
struct SavedVarConf
#	"A reference to the variable."
#	var :: VariableRef
	"Stores wether the variable was binary."
	was_bin :: Bool
	"Stores wether the variable was integer (not binary, nor continuous)."
	was_int :: Bool
	"Stores wether the variable was fixed."
	was_fixed :: Bool
	"If the variable was fixed, to which value they were fixed."
	fix_value :: Float64
	"Stores wether the variable had a lower bound."
	had_lb :: Bool
	"If the variable had a lower bound, the value of their lower bound."
	lb :: Float64
	"Stores wether the variable had an upper bound."
	had_ub :: Bool
	"If the variable had an upper bound, the value of their upper bound."
	ub :: Float64
end
export SavedVarConf

"""
    restore!(var :: VariableRef, c :: SavedVarConf) :: Nothing

If `var` type and/or bounds are different than the ones specified in `c`,
then change `var` type and/or bounds to adhere to `c`.
"""
function restore!(var :: VariableRef, c :: SavedVarConf) :: Nothing
	!c.was_bin && is_binary(var) && unset_binary(var)
	!c.was_int && is_integer(var) && unset_integer(var)
	c.was_bin && !is_binary(var) && set_binary(var)
	c.was_int && !is_integer(var) && set_integer(var)

	if c.was_fixed && (!is_fixed(var) || fix_value(var) != c.fix_value)
		fix(var, c.fix_value; force = true)
	else
		is_fixed(var) && unfix(var)
		if c.had_lb && (!has_lower_bound(var) || lower_bound(var) != c.lb)
			set_lower_bound(var, c.lb)
		end
		if c.had_ub && (!has_upper_bound(var) || upper_bound(var) != c.ub)
			set_upper_bound(var, c.ub)
		end
	end
	nothing
end
export restore!

"""
    SavedVarConf(var :: VariableRef) :: SavedVarConf

Creates a `SavedVarConf` struct from the configuration of the given variable.
Note that the VariableRef itself is not stored.
"""
function SavedVarConf(var :: VariableRef) :: SavedVarConf
	was_bin = is_binary(var)
	was_int = is_integer(var)
	was_fixed = is_fixed(var)
	had_lb = has_lower_bound(var)
	had_ub = has_upper_bound(var)
	return SavedVarConf(
		was_bin, was_int,
		was_fixed, (was_fixed ? fix_value(var) : 0.0),
		had_lb, (had_lb ? lower_bound(var) : 0.0),
		had_ub, (had_ub ? upper_bound(var) : 0.0)
	)
end

"""
    save_and_fix!(var, fix_value = 0.0) :: SavedVarConf

Fix the variable value and return its SavedVarConf before the fix.
"""
function save_and_fix!(var :: VariableRef, fix_value = 0.0)
	svc = SavedVarConf(var)
	fix(var, fix_value; force = true)
	return svc
end
export save_and_fix!

"""
    relax!(var) -> var

The `var` is made continuous. If the variable was binary, and had a lower
(upper) bound below zero (above one) it is replaced by zero (one).
"""
function relax!(var :: VariableRef)
	if is_integer(var)
		unset_integer(var)
	elseif is_binary(var)
		unset_binary(var)
		if !has_lower_bound(var) || lower_bound(var) < 0.0
			set_lower_bound(var, 0.0)
		end
		if !has_upper_bound(var) || upper_bound(var) < 0.0
			set_upper_bound(var, 1.0)
		end
	end

	return var
end
export relax!

"""
    save_and_relax!(var) :: SavedVarConf

Does the same as `relax!` but save the variable configuration first and
return it.

See also: [`relax!`](@ref)
"""
function save_and_relax!(var :: VariableRef)
	svc = SavedVarConf(var)
	if svc.was_int
		unset_integer(var)
	elseif svc.was_bin
		unset_binary(var)
		if !svc.had_lb || svc.lb < 0.0
			set_lower_bound(var, 0.0)
		end
		if !svc.had_ub || svc.ub > 1.0
			set_upper_bound(var, 1.0)
		end
	end

	return svc
end
export save_and_relax!

# Not used by the rest of Utilities, just made available to other modules.
include("Args.jl")

end # module
