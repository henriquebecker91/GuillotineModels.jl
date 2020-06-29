module CommandLine

# TODO: accept command-line arguments to set the instance parsing configuration

# External packages used.
import TimerOutputs.@timeit
using JuMP
import MathOptInterface
const MOI = MathOptInterface
import Dates
import Dates.@dateformat_str
import ArgParse
import ArgParse: ArgParseSettings, set_default_arg_group!, add_arg_group!
import ArgParse: @add_arg_table!, add_arg_table!

include("./SolversArgs.jl") # empty_configured_model, *parse_settings()
using .SolversArgs

import ..TIMER # The global module timer.
using ..InstanceReader
import ..build_model
using ..Utilities
using ..Utilities.Args
using ..PPG2KP, ..PPG2KP.Args
using ..Flow, ..Flow.Args
#using ..KnapsackPlates, ..KnapsackPlates.Args
import ..get_cut_pattern, ..to_pretty_str, ..simplify!
# Used inside read_build_solve_and_print in tandem with generic-time-limit.
import ..throw_if_timeout_now, ..TimeoutError

"""
    create_unprefixed_subset(prefix, p_args :: T) :: T

!!! **Internal use.**

Given some `prefix`, query the `accepted_arg_list(Val{Symbol(prefix)})`,
to know which arguments were prefixed this way, search for them (with
the prefix) in `p_args` and return a new `typeof(p_args)` object in
which there is only the searched key-value pairs but the keys are
changed to not have the prefix anymore.
"""
function create_unprefixed_subset(prefix, p_args :: T) :: T where {T}
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

!!! **Internal use.**

Given two integers and two integer arrays, uses `p_args` keys
`div-and-round-{nearest,up,down}` to either: (1) return them unmodified
if all keys have value one; (2) return a copy of them that is rounded
the specified way (no two keys may have a value different than one).
If a copy is returned, the type of the scalars and the element type of
the arrays is converted to typeof(L).
"""
@timeit TIMER function div_and_round_instance(L, W, l, w, p_args)
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
	L, W, l, w
end

# Read the instance, build the model, solve the model, and print related stats.
"""
    read_build_solve_and_print(pp)

!!! **Internal use.**

Given the parsed parameters (`pp`), read the instance file, build the model,
solve the model (unless `pp['do-not-solve']` is true), and print statistics
related to this process (unless `pp['no-csv-output']` is true).
The list of options recognized and implemented by this method is the
list returned by `generic_args()` (it also needs the required arguments
in `core_argparse_settings()`). The other arguments are solver or model
specific and are extracted and passed to their specific methods.
"""
@timeit TIMER function read_build_solve_and_print(pp) # pp stands for parsed parameters
	start :: Float64 = time() # seconds since epoch
	limit :: Float64 = pp["generic-time-limit"]
	instfname = pp["instfname"]
	!pp["no-csv-output"] && @show instfname

	N, L_, W_, l_, w_, p, d = InstanceReader.read_from_file(instfname)
	before_build_time = time()
	L, W, l, w = div_and_round_instance(L_, W_, l_, w_, pp)
	throw_if_timeout_now(start, limit)

	solver_pp = create_unprefixed_subset(pp["solver"], pp)
	m = empty_configured_model(Val(Symbol(pp["solver"])), solver_pp)
	throw_if_timeout_now(start, limit)

	model_pp = create_unprefixed_subset(pp["model"], pp)
	model_id = Val(Symbol(pp["model"]))
	bmr = build_model(
		model_id, m, d, p, l, w, L, W, model_pp
	)
	throw_if_timeout_now(start, limit)

	!isempty(pp["save-model"]) && !pp["no-csv-output"] &&
		@timeit TIMER "save_model" JuMP.write_to_file(m, pp["save-model"])
	throw_if_timeout_now(start, limit)

	p_ = max(L/minimum(l), W/minimum(w))
	num_vars = num_variables(m)
	num_constrs = num_all_constraints(m)
	if !pp["no-csv-output"]
		#@show time_to_build_model
		n = length(d)
		@show n
		n_ = sum(d)
		@show n_
		@show p_
		@show num_vars
		@show num_constrs
	end

	pp["do-not-solve"] && return nothing

	throw_if_timeout_now(start, limit)
	pp["no-csv-output"] || println("MARK_FINAL_GENERIC_SOLVE")
	output_name = "finished_model_solve"
	output_value = @elapsed begin
		@timeit TIMER output_name optimize_within_time_limit!(m, start, limit)
	end
	@assert termination_status(m) == MOI.OPTIMAL
	pp["no-csv-output"] || println("$output_name = $output_value")
	after_solve_time = time()
	build_and_solve_time = after_solve_time - before_build_time
	!pp["no-csv-output"] && @show build_and_solve_time
	#See comment above about TimerOutputs.
	#time_to_solve_model = TimerOutputs.time(get_defaulttimer(), "optimize!")

	# This needs to be done even if it is not printed to warm-start the JIT.
	obj_value = 0.0
	if primal_status(m) == MOI.FEASIBLE_POINT
		@timeit TIMER "stringfy_solutions" begin
			obj_value = objective_value(m)
			solution = get_cut_pattern(model_id, m, eltype(d), eltype(l), bmr)
			sol_str = "solution = $solution\n"
			sol_pretty_str = "PRETTY_STR_SOLUTION_BEGIN\n" * to_pretty_str(solution) *
				"\nPRETTY_STR_SOLUTION_END\n"
			sol_simple_pretty_sol = "SIMPLIFIED_PRETTY_STR_SOLUTION_BEGIN\n" *
				to_pretty_str(simplify!(deepcopy(solution))) *
				"\nSIMPLIFIED_PRETTY_STR_SOLUTION_END\n"
		end
	end
	if !pp["no-csv-output"]
		#@show time_to_solve_model
		if primal_status(m) == MOI.FEASIBLE_POINT
			@timeit TIMER "print_solutions" begin
				iob = IOBuffer()
				write(iob, sol_str)
				write(iob, sol_pretty_str)
				write(iob, sol_simple_pretty_sol)
				print(read(seekstart(iob), String))
			end
		end
		@show obj_value
		obj_bound = objective_bound(m)
		@show obj_bound
		stop_reason = termination_status(m)
		@show stop_reason
		stop_code = Int64(stop_reason)
		@show stop_code
	end

	return nothing
end

"""
    gen_prefixed_argparse_settings(solver_or_model_name) :: ArgParseSettings

!!! **Internal use.**

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
	ArgParse.add_arg_group!(
		s, "$(prefix)-specific Options", "$(prefix)-specific-options"
	)
	original_args = accepted_arg_list(Val(Symbol(prefix)))
	prefixed_args = map(original_args) do arg
		Arg(string(prefix) * "-" * arg.name, arg.default, arg.help)
	end
	for arg in prefixed_args
		ArgParse.add_arg_table!(s, arg)
	end
	set_default_arg_group!(s)
	s
end

"""
    core_argparse_settings() :: ArgParseSettings

!!! **Internal use.**

An ArgParseSettings with the three core positional arguments `model`,
`solver`, and `instfname`. They cannot be modeled as `Arg` objects because,
by design, all extra arguments must be options (i.e., be optional and preceded
by dashes).
"""
function core_argparse_settings() :: ArgParseSettings
	s = ArgParse.ArgParseSettings()
	ArgParse.add_arg_group!(s, "Core Parameters", "core_parameters")
	@add_arg_table! s begin
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
	set_default_arg_group!(s)
	s
end

"""
    generic_args() :: Vector{Arg}

!!! **Internal use.**

All the `Arg`s representing options that are independent from solver or model.
"""
function generic_args() :: Vector{Arg}
	return [
		Arg(
			"do-not-solve", false,
			"The model is built but not solved. A solver has yet to be specified. Note that just building a model may depend on using a solver over subproblems. Such uses of the solver are not disabled by this flag."
		),
		Arg(
			"generic-time-limit", 60.0 * 60.0 * 24.0 * 365.0,
			"Defines a time limit (in seconds) to be observed in the context of the model-agnostic proccess of reading, building, solving, and printing. Each model has to define its own flag to control the time inside the model building process (to be set independently, it is does not interact with this flag). `JuMP.set_time_limit_sec` is called over the model before starting to solve it, with the remaining time after reading instance and building the model (not the total time). If `--warm-jit` defines that there will be a mock run to warm the jit, the time limit applies to this mock run (i.e., if the mock run breaks the limit it an exception is raised) but the timer is reset after the mock, so the mock time is not counted for the 'real' run time limit. A `GuillotineModels.TimeoutError` is raised if the time limit is violated."
		),
		Arg(
			"save-model", "",
			"Save the model of the passed instance to the given filename (use the extension to define the format, allowed formats are described by the enumeration MathOptInterface.FileFormats.FileFormat)."
		),
		Arg(
			"warm-jit", "with-toy",
			"There are three valid options: 'no', 'with-toy', and 'yes'. If 'yes', the instance is solved two times, but the first time the arguments are changed to disable output (i.e., '--no-csv-output' and '--SOLVER-no-output' are passed, and if supported, '--MODEL-quiet' is passed too; in the GuillotineModels.TIMER the timings of this first run are under 'warm-jit'). If 'with-toy', it behaves similarly to 'yes' but instead of using the instance itself it uses a very small hardcoded instance (this helps a lot, but many procedures only called when the model has specific properties are not called, and therefore, not compiled). If 'no', the code is just run a single time."
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

!!! **Internal use.**

An ArgParseSettings containing all options (not positional arguments) that are
independent from the chosen solver or model.
"""
function generic_argparse_settings()
	s = ArgParse.ArgParseSettings()
	add_arg_group!(s, "Generic Options", "generic_options")
	for arg in generic_args()
		add_arg_table!(s, arg)
	end
	set_default_arg_group!(s)
	s
end

"""
    argparse_settings(models_list, solvers_list) :: ArgParseSettings

!!! **Internal use.**

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
	ArgParse.import_settings!(s, generic_argparse_settings())
	for sym in [solvers_list; models_list]
		ArgParse.import_settings!(s, gen_prefixed_argparse_settings(sym))
	end
	s
end

"""
    warn_if_changed_unused_values(p_args, models_list, solvers_list)

!!! **Internal use.**

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
			if unused_arg_idx === nothing
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

!!! **Internal use.**

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
	# The flag conflicts of the specific model and solver.
	Utilities.Args.throw_if_incompatible_options(
		Val(Symbol(p_args["model"])),
		create_unprefixed_subset(p_args["model"], p_args)
	)
	Utilities.Args.throw_if_incompatible_options(
		Val(Symbol(p_args["solver"])),
		create_unprefixed_subset(p_args["solver"], p_args)
	)
	throw_if_unrecognized(
		"warm-jit", p_args["warm-jit"], ("no", "with-toy", "yes")
	)
	option_name = "generic-time-limit"
	p_args[option_name] > 0.0 || throw(ArgumentError(
		"Option $(option_name) must be positive, but it is" *
		" $(p_args[option_name])."
	))
end

# Definition of the command line arguments.
"""
    parse_args(args, models_list, solvers_list) :: Dict{String, Any}

!!! **Internal use.**

Given a vector of the command-line arguments `args` and the lists of
available models and solvers, parse the arguments. If the `args`
refer to a model or solver not in `models_list` or `solvers_list`
exceptions may be thrown. If `args` just triggers the help message,
an empty `Dict` is returned.
"""
@timeit TIMER function parse_args(
	args, models_list, solvers_list
) :: Dict{String, Any}
	s = argparse_settings(models_list, solvers_list)
	p_args = ArgParse.parse_args(args, s)
	isnothing(p_args) && return Dict{String, Any}()
	warn_if_changed_unused_values(p_args, models_list, solvers_list)
	throw_if_incompatible_options(p_args)

	return p_args
end

"""
    toy_instance() :: String

Just return a string representing a small instance to be used by `warm_jit`
if `--warm-jit with-toy`.
"""
function toy_instance() :: String
	"""
	100	200
	2
	50     100    1     2
	50     200    1     1
	"""
end

"""
    warm_jit(p_args)

!!! **Internal use.**

Looks at `p_args["warm-jit"]` and implements what is said in the flag
description.

See also: [`toy_instance`](@ref), [`read_build_solve_and_print`](@ref)
"""
@timeit TIMER function warm_jit(p_args)
	p_args["warm-jit"] == "no" && return
	copyed_args = copy(p_args)
	# The temporary file is only needed if p_args["warm-jit"] == "with-toy",
	# however, if it is needed, then it needs to be valid for the rest of
	# this method scope.
	mktemp() do path, io
		if p_args["warm-jit"] == "with-toy"
				write(io, toy_instance())
				close(io)
				copyed_args["instfname"] = path
		end
		copyed_args["no-csv-output"] = true
		copyed_args[p_args["solver"] * "-" * "no-output"] = true
		valid_model_args = accepted_arg_list(Val(Symbol(p_args["model"])))
		if any(arg -> arg.name == "quiet", valid_model_args)
			copyed_args[p_args["model"] * "-" * "quiet"] = true
		end
		read_build_solve_and_print(copyed_args)
	end
	return
end

"""
    run(args = ARGS; implemented_models = [...], supported_solvers = [...])

Parse the command-line arguments and take appropriate action: prints help
or solves the problem instance using the specified solver, model, and
their options.

The parameters available (listed in the help message and actually parsed)
depend on the `implemented_models` and `supported_solvers`. The default values
are the models and solvers made available by the GuillotineModels package. If
you implement your own models or add support for more solvers you need to pass
them there for them to be considered (and if you want the old ones to keep
working you need to specify them also).

The best way to know everything this command is capable is to call:
```
import SUPPORTED_SOLVERS_YOU_HAVE_AVAILABLE
import GuillotineModels
GuillotineModels.CommandLine.run(
	["--help"];
	supported_solvers = [SUPPORTED_SOLVERS_YOU_HAVE_AVAILABLE]
)
```

Which will give you the help message. If you want to use it as an script
you just need to remove the `["--help"]` from the call.
"""
@timeit TIMER function run(
	args = ARGS;
	implemented_models :: Vector{Symbol} = [:Flow, :PPG2KP],#, :KnapsackPlates],
	supported_solvers :: Vector{Symbol} = [:CPLEX, :Gurobi, :Cbc, :GLPK]
)
	p_args = parse_args(args, implemented_models, supported_solvers)
	isempty(p_args) && return # Happens if called just for "--help".
	!p_args["no-csv-output"] && @show args
	!p_args["no-csv-output"] && @show p_args
	date_now = Dates.format(Dates.now(), dateformat"yyyy-mm-ddTHH:MM:SS")
	!p_args["no-csv-output"] && @show date_now
	warm_jit(p_args) # It checks `p_args["warm-jit"]` to check if it will run.
	# remaining code should not depend on this option
	delete!(p_args, "warm-jit")
	total_instance_time = @elapsed read_build_solve_and_print(p_args)
	!p_args["no-csv-output"] && @show total_instance_time

	return
end

end # module
