module Utilities

using DocStringExtensions # for TYPEDFIELDS

"""
    bits2idxs(bits, idx_type = Int) :: Vector{idx_type}

Given a `BitArray` or a `Vector{Bool}` return a vector of the indexes storing
true values.
"""
function bits2idxs(bits) :: Vector{Int}
	return bits2idxs(bits, Int) :: Vector{Int}
end
function bits2idxs(bits, ::Type{I}) :: Vector{I} where {I}
	n = count(bits)
	iszero(n) && return I[]
	idxs = Vector{I}(undef, n)
	i = 1
	for (j, v) in pairs(bits)
		if v
			idxs[i] = convert(I, j)
			i += 1
			i > n && break
		end
	end
	return idxs
end
export bits2idxs

"""
    gather_nonzero(vars, ::Type{D}, threshold = 1e-5, sol_idx = 1)

Given some valid JuMP.Model `vars`, return a list of all indexes in `vars`
in which the variable value rounded to nearest integer is non-zero, and
a list of the rounded values themselves.

The `::Type{D}` is the integer type for the rounded values. The `threshold`
parameter is used to give a warning if the difference between the extracted
value (Float64) and the rounded value is larger than it. The `sol_idx` is
passed to `JuMP.value(...; result = sol_idx)` to allow choosing which solution
of the model is to be queried.
"""
function gather_nonzero(
	vars, ::Type{D}, threshold = 1e-5, sol_idx = 1
) where {D}
	@assert threshold <= 0.5
	values = value.(vars; result = sol_idx)
	rounded_values = round.(D, values, RoundNearest)
	if !isone(threshold)
		diffs = abs.(values .- rounded_values)
		large_diffs = filter(>(threshold), diffs)
		if !isempty(large_diffs)
			@warn "The value of $(length(large_diffs)) model variables is" *
				" different from their nearest integer by at least $(threshold)," *
				" the largest difference is $(maximum(large_diffs)). This may mean:" *
				" the tolerances are too lax; the model was an LP model (in this" *
				" case set `threshold` to 0.5); there is some numerical error;" *
				" or some variables are not integers as they should."
		end
	end
	nz_bits = (!iszero).(rounded_values)
	idxs = bits2idxs(nz_bits)
	vals = deleteat!(rounded_values, .!nz_bits)
	return idxs, vals
end
export gather_nonzero

"""
    unify!(::Type{QT_TYPE}, a)

Apply `sort!` and `unique!` to array `a` and then returns a `Vector{QT_TYPE}`
with the corresponding quantity of value in `a` before compacting it.


```julia-repl
> a = [40, 10, 20, 10, 20, 30, 20];

> unify!(Int16, a)
4-element Array{Int16,1}:
 2
 3
 1
 1

> a
4-element Array{Int64,1}:
 10
 20
 30
 40

```
"""
function unify!(::Type{QT_TYPE}, a) where {QT_TYPE}
	isempty(a) && return QT_TYPE[]
	sort!(a)
	last_v = first(a)
	num_distinct_values = one(QT_TYPE)
	for v in a
		if v != last_v
			num_distinct_values += one(QT_TYPE)
			last_v = v
		end
	end
	num_distinct_values == length(a) && return ones(QT_TYPE, num_distinct_values)
	qts = zeros(QT_TYPE, num_distinct_values)
	qts_idx = 1
	last_v = first(a)
	for v in a
		if v != last_v
			last_v = v
			qts_idx += 1
		end
		qts[qts_idx] += one(QT_TYPE)
	end

	unique!(a)
	return qts
end
export unify!

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

function flush_all_output()
	flush(stdout)
	flush(stderr)
end
export flush_all_output

using JuMP

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
#=
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
=#

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

# Internal structure. It does not validate its input but it should be clear
# that the values of 'delete_X/unfix' and 'set_X/fix' cannot be both true (both
# may be false, this means that nothing is done). `fix` and any `set_` cannot
# be both true either. The idea is to store knowledge of "what to change"
# that depends on both the SavedVarConf and the current state of the variable.
# With this structure, all variables may be queried and compared to the saved
# configurations before the first variable is changed. See the comment at
# `restore!`.
struct _VarConfDiff
	unset_bin :: Bool
	set_bin :: Bool
	unset_int :: Bool
	set_int :: Bool
	unfix :: Bool
	fix :: Bool
	fix_value :: Float64
	delete_lb :: Bool
	set_lb :: Bool
	lb_value :: Float64
	delete_ub :: Bool
	set_ub :: Bool
	ub_value :: Float64
end

function _diff(var :: VariableRef, c :: SavedVarConf) :: _VarConfDiff
	is_bin = is_binary(var)
	is_int = (is_bin ? false : is_integer(var))
	has_fix = is_fixed(var)
	has_lb = (has_fix ? false : has_lower_bound(var))
	has_ub = (has_fix ? false : has_upper_bound(var))
	return _VarConfDiff(
		!c.was_bin & is_bin      # unset_bin
		, c.was_bin & !is_bin    # set_bin
		, !c.was_int & is_int    # unset_int
		, c.was_int & !is_int    # set_int
		, !c.was_fixed & has_fix # unfix
		, c.was_fixed && (!has_fix || fix_value(var) != c.fix_value) # fix
		, c.fix_value # fix_value
		, !c.had_lb && has_lb # delete_lb
		, c.had_lb && (!has_lb || lower_bound(var) != c.lb) # set_lb
		, c.lb # lb_value
		, !c.had_ub && has_ub # delete_ub
		, c.had_ub && (!has_ub || upper_bound(var) != c.ub) # set_ub
		, c.ub # ub_value
	)
end

function _apply(var :: VariableRef, diff :: _VarConfDiff) :: Nothing
	# Unset properties before setting ones that cannot coexist with the first.
	diff.unset_bin && unset_binary(var)
	diff.unset_int && unset_integer(var)
	diff.set_bin && set_binary(var)
	diff.set_int && set_integer(var)
	diff.delete_lb && delete_lower_bound(var)
	diff.delete_ub && delete_upper_bound(var)
	if diff.fix
		# Alternatively, we could use `force = true`, but we already have the
		# deletion information queried, and deleted bounds if they existed above
		# so we can call fix without `force = true`.
		fix(var, diff.fix_value)
	else
		diff.unfix && unfix(var)
		diff.set_lb && set_lower_bound(var, diff.lb_value)
		diff.set_ub && set_upper_bound(var, diff.ub_value)
	end
	return
end

# Because of how Gurobi works this method need be be written carefully
# (see https://github.com/JuliaOpt/Gurobi.jl/pull/301). Basically, we
# want to query every attribute from the model before setting any of them.
# Also, broadcasting this method is not a good idea, a batch version is
# necessary.
"""
    restore!(var :: VariableRef, c :: SavedVarConf) :: Nothing
    restore!(vars :: Vector{...}, cs :: Vector{...}) :: Nothing

If `var` type and/or bounds are different than the ones specified in `c`,
then change `var` type and/or bounds to adhere to `c`.
"""
function restore!(var :: VariableRef, c :: SavedVarConf) :: Nothing
	_apply(var, _diff(var, c))
	return
end
export restore!

function restore!(
	vars :: AbstractVector,
	cs :: AbstractVector
) :: Nothing
	@assert length(vars) == length(cs)
	# do not nest the broadcasts, we do not want them to fuse
	ds = _diff.(vars, cs)
	_apply.(vars, ds)
	return
end

"""
    SavedVarConf(var :: VariableRef) :: SavedVarConf

Creates a `SavedVarConf` struct from the configuration of the given variable.
Note that the VariableRef itself is not stored.
"""
function SavedVarConf(var :: VariableRef) :: SavedVarConf
	was_bin = is_binary(var)
	was_int = (was_bin ? false : is_integer(var))
	was_fixed = is_fixed(var)
	had_lb = (was_fixed ? false : has_lower_bound(var))
	had_ub = (was_fixed ? false : has_upper_bound(var))
	return SavedVarConf(
		was_bin, was_int,
		was_fixed, (was_fixed ? fix_value(var) : 0.0),
		had_lb, (had_lb ? lower_bound(var) : 0.0),
		had_ub, (had_ub ? upper_bound(var) : 0.0)
	)
end

function _relax!(var :: VariableRef, svc :: SavedVarConf) :: Nothing
	if svc.was_int
		unset_integer(var)
	elseif svc.was_bin
		unset_binary(var)
		(!svc.had_lb || svc.lb < 0.0) && set_lower_bound(var, 0.0)
		(!svc.had_ub || svc.ub > 1.0) && set_upper_bound(var, 1.0)
	end
	return
end

"""
    relax!(var) :: SavedVarConf
    relax!(vars) :: Vector{SavedVarConf}

The `var` is made continuous. If the variable was binary, and had a lower
(upper) bound below zero (above one) it is replaced by zero (one).
"""
function relax!(var :: VariableRef) :: SavedVarConf
	svc = SavedVarConf(var)
	_relax!(var, svc)
	return svc
end
export relax!

function relax!(vars :: AbstractVector) :: Vector{SavedVarConf}
	# do not nest the broadcasts, we do not want them to fuse
	svcs = SavedVarConf.(vars)
	_relax!.(vars, svcs)
	return svcs
end

# Not used by the rest of Utilities, just made available to other modules.
include("Args.jl")

using Statistics
function vector_summary(v :: AbstractVector)
	println(summary(v))
	if !isempty(v)
		q0_to_100 = quantile(v, 0.0:0.1:1.0)
		@show q0_to_100
		@show mean(v)
		@show std(v)
		@show var(v)
		num_zeros = sum(isapprox.(0.0, v; atol = 1e-6))
		@show num_zeros
	end
end
export vector_summary

end # module
