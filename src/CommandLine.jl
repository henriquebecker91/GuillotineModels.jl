module CommandLine

# TODO: accept command-line arguments to set the instance parsing configuration

# External packages used.
using TimerOutputs
using JuMP
import ArgParse
import ArgParse: ArgParseSettings, set_default_arg_group, add_arg_group
import ArgParse: @add_arg_table, add_arg_table

include("./SolversArgs.jl") # empty_configured_model, *parse_settings()
using .SolversArgs

using ..InstanceReader
import ..build_model
using ..Utilities
using ..Utilities.Args
using ..PPG2KP, ..PPG2KP.Args
using ..Flow, ..Flow.Args
using ..KnapsackPlates, ..KnapsackPlates.Args

"""
    create_unprefixed_subset(prefix, p_args :: T) :: T

Given some `prefix`, query the `accepted_arg_list(Val{Symbol(prefix)})`,
to know which arguments were prefixed this way, search for them (with
the prefix) in `p_args` and return a new `typeof(p_args)` object in
which there is only the searched key-value pairs but the keys are
changed to not have the prefix anymore.
"""
function create_unprefixed_subset(prefix, p_args :: T) :: T
	subset = empty(p_args)
	prefix_id = Val(Symbol(prefix))
	prefix_options_names = getfield.(accepted_arg_list(prefix_id), :name)
	for name in prefix_options_names
		prefixed_name = prefix * "-" * name
		@assert haskey(p_args, prefixed_name)
		subset[name] = p_args[prefixed_name]
	end
	subset
end

"""
    div_and_round_instance(L, W, l, w, p_args) -> L', W', l', w'

Given two integers and two integer arrays, uses `p_args` keys
`div-and-round-{nearest,up,down}` to either: (1) return them unmodified
if all keys have value one; (2) return a copy of them that is rounded
the specified way (no two keys may have a value different than one).
If a copy is returned, the type of the scalars and the element type of
the arrays is converted to typeof(L).
"""
function div_and_round_instance(L, W, l, w, p_args)
	@timeit "div_and_round_instance" begin
	# assert explanation: at least two of the three flags are disabled (i.e., 
	# have value one)
	@assert sum(isone.((
		p_args["div-and-round-nearest"],
		p_args["div-and-round-up"],
		p_args["div-and-round-down"]
	))) >= 2

	if p_args["div-and-round-nearest"] != 1
		factor = p_args["div-and-round-nearest"]
		roundmode = RoundNearest
	elseif p_args["div-and-round-up"] != 1
		factor = p_args["div-and-round-up"]
		roundmode = RoundUp
	elseif p_args["div-and-round-down"] != 1
		factor = p_args["div-and-round-down"]
		roundmode = RoundDown
	else
		factor = 1
	end
	if factor != 1
		S = eltype(L)
		L = convert(S, round(L / factor, roundmode))
		W = convert(S, round(W / factor, roundmode))
		l = convert.(S, round.(l ./ factor, roundmode))
		w = convert.(S, round.(w ./ factor, roundmode))
	end
	end # timeit
	L, W, l, w
end

# Read the instance, build the model, solve the model, and print related stats.
"""
    read_build_solve_and_print(pp)

Given the parsed parameters (`pp`), read the instance file, build the model,
solve the model (unless `pp['do-not-solve']` is true), and print statistics
related to this process (unless `pp['no-csv-output']` is true).
The list of options recognized and implemented by this method is the
list returned by `generic_args()` (it also needs the required arguments
in `core_argparse_settings()`). The other arguments are solver or model
specific and are extracted and passed to their specific methods.
"""
function read_build_solve_and_print(pp) # pp stands for parsed parameters
	if !pp["no-csv-output"]
		println("instfname = $(pp["instfname"])")
	end

	N, L_, W_, l_, w_, p, d = InstanceReader.read_from_file(pp["instfname"])
	L, W, l, w = div_and_round_instance(L_, W_, l_, w_, pp)

	solver_pp = create_unprefixed_subset(pp["solver"], pp)
	@timeit "empty_configured_model" begin
	m = empty_configured_model(Val(Symbol(pp["solver"])), solver_pp)
	end # timeit

	@timeit "build_model" begin
	model_pp = create_unprefixed_subset(pp["model"], pp)
	build_model_return = build_model(
		Val(Symbol(pp["model"])), m, d, p, l, w, L, W, model_pp
	)
	end # @timeit
	#@show build_model_return
	#= This does not work unless we give the full path, what is a load of shit.
	# Needs to be fixed either by: (1) PR to the package solving the problem;
	# (2) using an @elapsed inside the block (or outside it); (3) hacking the
	# package internals to know the already existing stack of sections.
	time_to_build_model = TimerOutputs.time(
		TimerOutputs.get_defaulttimer()["build_model"]
	)
	=#

	pp["save-model"] && !pp["no-csv-output"] &&
		@timeit "save_model" save_model(m, "./$(basename(instfname)).mps")

	if !pp["no-csv-output"]
		#@show time_to_build_model
		n = length(d)
		@show n
		n_ = sum(d)
		@show n_
		p_ = max(L/minimum(l), W/minimum(w))
		@show p_
		num_vars = num_variables(m)
		@show num_vars
		num_constrs = num_all_constraints(m)
		@show num_constrs
		#= TODO: make specific prints for each model inside their module and
		# call them there?
		if !pp["flow-model"] && !pp["break-hvcut-symmetry"]
			num_plates = length(pli_lwb)
			@show num_plates
			ps = m[:picuts]
			num_picuts = length(ps)
			@show num_picuts
			cm = m[:cuts_made]
			num_cuts_made = length(cm)
			@show num_cuts_made
		end
		=#
	end

	pp["do-not-solve"] && return nothing

	@timeit "optimize!" optimize!(m)
	#See comment above about TimerOutputs.
	#time_to_solve_model = TimerOutputs.time(get_defaulttimer(), "optimize!")

	if !pp["no-csv-output"]
		#@show time_to_solve_model
		if primal_status(m) == MOI.FEASIBLE_POINT
			obj_value = objective_value(m)
		else
			obj_value = 0
		end
		@show obj_value
		if !pp["relax2lp"]
			obj_bound = objective_bound(m)
			@show obj_bound
		end
		stop_reason = termination_status(m)
		@show stop_reason
		stop_code = Int64(stop_reason)
		@show stop_code
	end

	return nothing
end

"""
    gen_prefixed_argparse_settings(solver_or_model_name) :: ArgParseSettings

Builds and returns an ArgParseSettings object representing all options
of the given solver or model already prefixed with its name.

The `solver_or_model_name` may be a Symbol or String (none is prefered),
and the method works because the solver or model implements
`accepted_arg_list(Val{Symbol(solver_or_model_name)})`
"""
function gen_prefixed_argparse_settings(
	solver_or_model_name :: Union{String, Symbol}
) :: ArgParse.ArgParseSettings
	prefix = solver_or_model_name # long name only for documentation
	s = ArgParseSettings()
	ArgParse.add_arg_group(
		s, "$(prefix)-specific Options", "$(prefix)-specific-options"
	)
	original_args = accepted_arg_list(Val(Symbol(prefix)))
	prefixed_args = map(original_args) do arg
		Arg(string(prefix) * "-" * arg.name, arg.default, arg.help)
	end
	for arg in prefixed_args
		ArgParse.add_arg_table(s, arg)
	end
	set_default_arg_group(s)
	s
end

"""
    core_argparse_settings() :: ArgParseSettings

An ArgParseSettings with the three core positional arguments `model`,
`solver`, and `instfname`. They cannot be modeled as `Arg` objects by design
all extra arguments must be options (i.e., be optional and preceded by dashes).
"""
function core_argparse_settings() :: ArgParseSettings
	s = ArgParse.ArgParseSettings()
	ArgParse.add_arg_group(s, "Core Parameters", "core_parameters")
	@add_arg_table s begin
		"model"
			help = "Model or solution procedure to be used (case sensitive, ex.: Flow, KnapsackPlates, PPG2KP). Required."
			arg_type = String
		"solver"
			help = "Solver to be used if necessary (case-sensitive, ex.: Cbc, CPLEX, Gurobi). Required, even if --do-not-solve is specified."
			arg_type = String
		"instfname"
			help = "The path to the instance to be solved."
			arg_type = String
	end
	set_default_arg_group(s)
	s
end

"""
    generic_args() :: Vector{Arg}

All the `Arg`s representing options that are independent from solver or model.
"""
function generic_args() :: Vector{Arg}
	return [
		Arg(
			"do-not-solve", false,
			"The model is build but not solved. A solver has yet to be specified. Note that just building a model may depend on using a solver over subproblems. Such uses of the solver are not disabled by this flag."
		),
		Arg( # TODO: check if it can be implemented generically
			"save-model", false,
			"Save the model of each problem instance to the working directory ('./instance_name.mps'). Uses MPS format."
		),
		Arg(
			"do-not-mock-for-compilation", false,
			"To avoid counting the compilation time, a small hardcoded mock instance is solved before (to warm the JIT), this flag disables this mock solve."
		),
		Arg(
			"no-csv-output", false,
			"Disable some extra output that the author of this package uses to assemble CSVs and/or debug. Also disables --save-model."
		),
		Arg(
			"relax2lp", false,
			"Integer and binary variables become continuous."
		),
		Arg(
			"div-and-round-nearest", 1,
			"Divide the instances lenght and width (but not profit) by the passed factor and round them to nearest (the model answer becomes a GUESS, not a valid primal heuristic, nor a valid bound)."
		),
		Arg(
			"div-and-round-up", 1,
			"Divide the instances lenght and width (but not profit) by the passed factor and round them up (the model becomes a PRIMAL HEURISTIC)."
		),
		Arg(
			"div-and-round-down", 1,
			"Divide the instances lenght and width (but not profit) by the passed factor and round them down (the model becomes an OPTIMISTIC GUESS, A VALID BOUND)."
		)
	]
end

"""
    generic_argparse_settings() :: ArgParseSettings

An ArgParseSettings containing all options (not positional arguments) that are
independent from the chosen solver or model.
"""
function generic_argparse_settings()
	s = ArgParse.ArgParseSettings()
	add_arg_group(s, "Generic Options", "generic_options")
	for arg in generic_args()
		add_arg_table(s, arg)
	end
	set_default_arg_group(s)
	s
end

"""
    argparse_settings(models_list, solvers_list) :: ArgParseSettings

Create an ArgParseSettings which includes the core arguments, generic
arguments, and all arguments from models and solvers available (prefixed by
their model or solver name).

Options of all models and solvers are added even if just one model is selected
at a time. The reasons for this are: (1) we do not know the model before the
parsing unless we manipulate ARGS directly; (2) if the unused solvers and
models are not included in the ArgParseSettings they do not appear in the help
message.
"""
function argparse_settings(
	models_list :: Vector{Symbol}, solvers_list :: Vector{Symbol}
) :: ArgParseSettings
	s = core_argparse_settings()
	ArgParse.import_settings(s, generic_argparse_settings())
	for sym in [solvers_list; models_list]
		ArgParse.import_settings(s, gen_prefixed_argparse_settings(sym))
	end
	s
end

"""
    warn_if_changed_unused_values(p_args, models_list, solvers_list)

Gives warning messages if `p_args` has a value different from the default
for an option of a model or solver that is not the one used. In other words,
help a distracted user to not keep thinking it has passed an option to the
used solver/model when they have not. Many solvers have the same option
but with a different prefix, it is easy to change the solver used and forget
to also change the prefix in the parameter options.

It is actually impossible to know if an argument was passed or not by the
command-line after they are parsed by ArgParse because all of the arguments
must have a default value and if the argument is not present in the
command-line their default value is placed in the dict by ArgParse. This is
the reason we just check if non-used options have values different from
default, instead of just checking if they were provided.

This should work out-of-the-box for any third-party models or solvers
given they implement their own version of `accepted_arg_list` and also
have their identifying symbol passed in either `models_list` or
`solvers_list`.
"""
function warn_if_changed_unused_values(p_args, models_list, solvers_list)
	core_arg_names = ["solver", "model", "instfname"]
	generic_arg_names = getfield.(generic_args(), :name)

	model = p_args["model"]
	model_id = Val(Symbol(model))
	model_arg_names = getfield.(accepted_arg_list(model_id), :name)
	model_arg_names = map(name -> model * "-" * name, model_arg_names)
	# the same as the block above just with "model" replaced by "solver"
	solver = p_args["solver"]
	solver_id = Val(Symbol(solver))
	solver_arg_names = getfield.(accepted_arg_list(solver_id), :name)
	solver_arg_names = map(name -> solver * "-" * name, solver_arg_names)

	expected_arg_names = [
		model_arg_names; core_arg_names; generic_arg_names; solver_arg_names
	]
	sort!(expected_arg_names)

	used_groups = Symbol.([p_args["model"], p_args["solver"]])
	unused_groups = setdiff([models_list; solvers_list], used_groups)
	unused_prefixed_args = Arg[]
	for unused_group in unused_groups
		args = accepted_arg_list(Val(unused_group))
		prefixed_args = map(args) do arg
			Arg(string(unused_group) * "-" * arg.name, arg.default, arg.help)
		end
		append!(unused_prefixed_args, prefixed_args)
	end
	unused_prefixed_args = collect(Iterators.flatten(unused_prefixed_args))

	for arg_name in keys(p_args)
		if isempty(searchsorted(expected_arg_names, arg_name))
			unused_arg_idx = findfirst(a -> a.name == arg_name, unused_prefixed_args)
			if isnothing(unused_arg_idx)
				@error(
					"ArgParse accepted option --$(arg_name) but it is not a" *
					" recognized option. This should not be possible."
				)
			else
				unused_arg = unused_prefixed_args[unused_arg_idx]
				unused_arg.default == p_args[arg_name] && continue
				@warn(
					"Option --$(unused_arg.name) has a different value than the" *
					" default, but it is not from the solver or model selected." *
					" You probably did use an option of an unselected solver or model."
				)
			end
		end
	end

	return nothing
end

"""
    throw_if_incompatible_options(p_args)

Check the already parsed arguments and test if options that are incompatible
with each other were used.
"""
function throw_if_incompatible_options(p_args)
	# Generic Flags conflicts
	num_rounds = sum(.! isone.(map(flag -> p_args[flag],
		["div-and-round-nearest", "div-and-round-up", "div-and-round-down"]
	)))
	num_rounds > 1 && @error(
		"only one of --div-and-round-{nearest,up,down} may be passed at the" *
		" same time (what the fuck you expected?)"
	)
	# The check below needs to be here because uses a generic flag and a specific
	# flag at same time.
	p_args["PPG2KP-final-pricing"] && p_args["relax2lp"] && @error(
		"The flags --final-pricing and --relax2lp should not be used together;" *
		" it is not clear what they should do, and the best interpretation" *
		" (solving the relaxed model and doing the final pricing, without" *
		" solving the unrelaxed reduced model after) is not specially useful" *
		" and need extra code to work that is not worth it."
	)
	Utilities.Args.throw_if_incompatible_options(
		Val(Symbol(p_args["model"])),
		create_unprefixed_subset(p_args["model"], p_args)
	)
	Utilities.Args.throw_if_incompatible_options(
		Val(Symbol(p_args["solver"])),
		create_unprefixed_subset(p_args["solver"], p_args)
	)
end

# Definition of the command line arguments.
"""
    parse_args(args, models_list, solvers_list) :: Dict{String, Any}

Given a vector of the command-line arguments `args` and the lists of
available models and solvers, parse the arguments. If the `args`
refer to a model or solver not in `models_list` or `solvers_list`
exceptions may be thrown.
"""
function parse_args(args, models_list, solvers_list) :: Dict{String, Any}
	@timeit "parse_args" begin
	s = argparse_settings(models_list, solvers_list)
	p_args = ArgParse.parse_args(args, s) :: Dict{String, Any}
	warn_if_changed_unused_values(p_args, models_list, solvers_list)
	throw_if_incompatible_options(p_args)
	end # @timeit

	return p_args
end

"""
    mock_instance() :: String

Just return a string of a small instance to be used by `mock_for_compilation`
if `--do-not-mock-for-compilation` is not specified.
"""
function mock_instance() :: String
	"""
	100	200
	2
	50     100    1     2
	50     200    1     1
	"""
end

"""
    mock_for_compilation(p_args)

Creates a temporary file with the content returned by `mock_instance` and pass
to `read_build_solve_and_print` a copy of `p_args` pointing to this instance
with `no-csv-output` and `USED_SOLVER_NAME-no-output` as true.

See also: [`mock_instance`](@ref), [`read_build_solve_and_print`](@ref)
"""
function mock_for_compilation(p_args)
	@timeit "mock_for_compilation" begin
	mktemp() do path, io
		copyed_args = copy(p_args)
		write(io, mock_instance())
		close(io) # it will be opened again inside read_build_solve_and_print
		copyed_args["instfname"] = path
		copyed_args["no-csv-output"] = true
		copyed_args[p_args["solver"] * "-" * "no-output"] = true
		read_build_solve_and_print(copyed_args)
	end
	end # @timeit
end

# Parse the command line arguments, and call the solve for each instance.
function run(
	args = ARGS;
	implemented_models = [:Flow, :PPG2KP],#, :KnapsackPlates],
	supported_solvers = [:CPLEX, :Gurobi, :Cbc, :GLPK]
)
	@timeit "run" begin
	#@show args
	p_args = parse_args(args, implemented_models, supported_solvers)
	#@show p_args
	!p_args["do-not-mock-for-compilation"] && mock_for_compilation(p_args)
	# remaining code should not depend on this option
	delete!(p_args, "do-not-mock-for-compilation")
	total_instance_time = @elapsed read_build_solve_and_print(p_args)
	@show total_instance_time
	end # timeit
end

end # module
