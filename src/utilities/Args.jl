"""
Submodule of Utilities that groups the utilities related to command-line
parsing.
"""
module Args

using DocStringExtensions # for TYPEDFIELDS
using ArgParse
export Arg, accepted_arg_list, throw_if_incompatible_options
export create_normalized_arg_subset

"""
An argument with `name`, `default` value, and `help` message.

$(TYPEDFIELDS)

"""
struct Arg{T}
  "The name with no initial double dashes but with dashes instead of spaces."
	name :: String
	"The default value if the argument is not specified. Use `false` for flags."
	default :: T
	"The help message to be displayed about this argument."
	help :: String
end

# The re-definition of Base.iterate is necessary to allow the flattening
# of arrays of arrays of Args.
Base.iterate(a :: Arg) = (a, nothing)
Base.iterate(a :: Arg, nothing) = nothing

"""
    accepted_arg_list(::Val{T}) :: Vector{Arg} where {T}

Generic error fallback. Any model or solver to be supported by
`GuillotineModules.SolversArgs.run` should implement their own version of this
method (replacing the `T` in `::Val{T}` by a Symbol identifying the solver
package or the name of the model). Look at module SolversArgs source for
examples of the implementation for a solver, and at PPG2KP module source for
examples for a model.

Gives the list of accepted options by some model or solver. The options have
either a boolean default (in this case they take no argument, the presence of
the option just flips the default value) or a non-boolean default (in this
case, if they are passed in the command-line, they must have a parameter; if
they are not passed, the default value is used).  Every solver implementation
must support the `no-output` and preferably the `raw-argument` option too. The
name of the supported solvers and implemented methods need to be passed to the
`run` method for them to be considered by it.
"""
function accepted_arg_list(::Val{T}) :: Vector{Arg} where {T}
	@error(
		"Solver " * string(T) * " is not supported (i.e., there is not" *
		" an implementation of Julia method `GuillotineModels.Utilities." *
		".Args.accepted_arg_list(::Val{Symbol(\"$(string(T))\")})` for it) or" *
		" the solver package was not imported before this method was called."
	)
end

"""
    throw_if_incompatible_options(::Val{T}, p_args) where {T}

Generic error fallback. Any model or solver to be supported by
`GuillotineModels.SolversArgs.run` should implement their own version of this
method (replacing the `T` in `::Val{T}` by a Symbol identifying the solver
package or the name of the model). Look at module SolversArgs source for
examples of the implementation for a solver, and at PPG2KP module source for
examples for a model. Every solver implementation must support the `no-output`
option, and there is no need to guarantee that conflicts with involving
the `raw-argument` option are detected.

Often is specialized to an empty method, as is not so often that solvers
or models have conflicting options.
"""
function throw_if_incompatible_options(
	::Val{T}, p_args
) where {T}
	@error (
		"A specialized method for $(T) should exist, but instead this " *
		" generic error fallback was called."
	)
end

"""
ArgParse.add_arg_table!(settings :: ArgParseSettings, arg :: Arg)

A specialization of `ArgParse.add_arg_table!` to transform `Arg` objects into
options of `ArgParseSettings`. Boolean arguments become `:store_{true|false}`
options (depending on the default) and non-boolean arguments become
`:store_arg` options with a default (and enforcing the same arg_type of the
default).
"""
function ArgParse.add_arg_table!(settings :: ArgParseSettings, arg :: Arg)
	if isa(arg.default, Bool)
		conf = Dict{Symbol, Any}(
			:help => arg.help,
			# default == true means that passing it as a flag turns it false
			:action => (arg.default ? :store_false : :store_true)
		)
	else
		conf = Dict{Symbol, Any}(
			:help => arg.help,
			:arg_type => typeof(arg.default),
			:default => arg.default,
			:action => :store_arg
		)
	end
	# only long options are accepted, enforced clarity
	ArgParse.add_arg_table!(settings, "--" * arg.name, conf)
end

"""
    create_normalized_arg_subset(p_args, selected :: Vector{Arg})

Create a new `Dict` with keys equal to `selected` names and values equal to
`selected` defaults, except that if `p_args` has a key with the same name
as an argument then the value in `p_args` is used instead the default.

It shows an warning for every key in `p_args` that is not a name in
`selected`.
"""
function create_normalized_arg_subset(
	p_args, selected :: Vector{Arg}
)
	new_dict = empty(p_args)
	for arg in selected
		if haskey(p_args, arg.name)
			new_dict[arg.name] = p_args[arg.name]
		else
			new_dict[arg.name] = arg.default
		end
	end
	sel_keys = sort(getfield.(selected, :name))
	old_keys = keys(p_args)
	for k in old_keys
		isempty(searchsorted(sel_keys, k)) && @warn(
			"Key $(k) is not recognized, it will not be used."
		)
	end

	return new_dict
end
end # submodule Args
