module CommandLine

# TODO: accept command-line arguments to set the instance parsing configuration

# External packages used.
using UnPack: @unpack
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
using ..Data
import ..build_model
import ..CutPattern, ..get_cut_pattern
import ..BuildStopReason, ..BUILT_MODEL, ..FOUND_OPTIMUM
using ..Utilities
using ..Utilities.Args
using ..PPG2KP, ..PPG2KP.Args
using ..Flow, ..Flow.Args
#using ..KnapsackPlates, ..KnapsackPlates.Args
import ..solution_value, ..to_pretty_str, ..simplify!
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
    round_instance(instance, factor, roundmode)

Return a new instance with all size-related fields multiplied by `factor`.

The `roundmode` (a `Base.RoundingMode`) specifies how the instances should
be rounded back to their original types (as the `factor` is expected to
be a `Float64` and not an `Integer`).

This function should be extended for new instance types; if you
want to use `--round-{nearest,up,down}` flags of
`GuillotineModels.CommandLine.run` for new instance types just create a
method that takes an instance of the new type, a `Float64`, and a
`Base.RoundingMode`).

This function is called with this signature only if at least one option
among `round-{nearest,up,down}` has a value different than one.

!!! Unmodified arrays (like `d`) may be shared between old and new instances.
"""
function round_instance(
	instance :: G2KP{D, S, P}, factor, roundmode
) where {D, S, P}
	@unpack L, W, l, w, d, p = instance
	L_, W_, l_, w_ = _round_instance(L, W, l, w, factor, roundmode)
	return G2KP{D, S, P}(L_, W_, l_, w_, d, p)
end

function round_instance(
	instance :: SLOPP{D, S, P}, factor, roundmode
) where {D, S, P}
	@unpack L, W, l, w, dlb, dub, p = instance
	L_, W_, l_, w_ = _round_instance(L, W, l, w, factor, roundmode)
	return SLOPP{D, S, P}(L_, W_, l_, w_, dlb, dub, p)
end

function _round_instance(
	L :: S, W :: S, l :: Vector{S}, w :: Vector{S}, factor, roundmode
) where {S}
	L_ = convert(S, round(L * factor, roundmode)) :: S
	W_ = convert(S, round(W * factor, roundmode)) :: S
	l_ = convert.(S, round.(l .* factor, roundmode)) :: Vector{S}
	w_ = convert.(S, round.(w .* factor, roundmode)) :: Vector{S}
	return L_, W_, l_, w_
end

"""
    round_instance(instance, p_args) -> old_or_new_instance

!!! **Internal use.**

Given a recognized instance type, uses `p_args` keys
`round-{nearest,up,down}` to either: (1) return them unmodified
if all keys have value one; (2) return a copy of them that is multiplied
by the ratio and rounded the specified way (no two keys may have a value
different than one).
"""
@timeit TIMER function round_instance(instance, p_args)
	# assert explanation: at least two of the three flags are disabled (i.e.,
	# have value one)
	@assert sum(isone.((
		p_args["round-nearest"],
		p_args["round-up"],
		p_args["round-down"]
	))) >= 2

	if p_args["round-nearest"] != 1
		factor = p_args["round-nearest"] :: Float64
		roundmode = RoundNearest
	elseif p_args["round-up"] != 1
		factor = p_args["round-up"] :: Float64
		roundmode = RoundUp
	elseif p_args["round-down"] != 1
		factor = p_args["round-down"] :: Float64
		roundmode = RoundDown
	else
		factor = 1
	end
	if factor != 1
		instance = round_instance(
			instance, factor, roundmode
		) :: typeof(instance)
	end
	return instance
end

# Fallback for unknown instance types.
function print_instance_stats(instance, verbose :: Bool)
	println("Start of the fallback output of `print_instance_stats`.")
	show(IOContext(stdout, :limit => false), instance)
	println("\nEnd of the fallback output of `print_instance_stats`.")
	return
end

function print_instance_stats(
	instance :: G2KP{D, S, P}, verbose :: Bool
) where {D, S, P}
	@unpack L, W, l, w, d = instance
	p_ = max(L/minimum(l), W/minimum(w))
	n_ = sum(d)
	n = length(d)
	if verbose
		@show L
		@show W
		@show n_
		@show n
		println("min_l = $(minimum(l))")
		println("min_w = $(minimum(w))")
		@show p_
		show(IOContext(stdout, :limit => false), instance)
		println()
	end
	return
end

# Get every $<something> in which the $ is not preceded by @
const TEMPLATE_REGEX = r"(?<!@)\$<([^>]+)>"
const ESCAPED_REGEX = r"@(\$<[^>]+>)"

# Treats `instance_file` specially.
function expand_template(
	template :: AbstractString, dict :: Dict{<:AbstractString, <:Any}
)
	expanded = replace(template, TEMPLATE_REGEX => function(s)
		# The below is byte indexing, not character indexing, but it is safe for
		# UFT-8 strings because '$<' are two bytes and '>' is another byte.
		key = s[3:end-1]
		if key == "instance_file"
			basename(dict["instance_path"])
		elseif haskey(dict, key)
			value = dict[key]
			if isa(value, Bool)
				string(convert(Int, value))
			else
				string(value)
			end
		else
			error("The template '\$<$key>' is invalid. There is no such parameter.")
		end
	end)
	cleaned = replace(expanded, ESCAPED_REGEX => s"\1")

	return cleaned
end

# Serves both as a function barrier and the step to be executed for each
# instance in the instance(s) file.
function build_solve_and_print(problem, formulation, instance_, pp, timings)
	start = time() :: Float64 # seconds since epoch
	limit = pp["generic-time-limit"] :: Float64
	verbose = !(pp["no-csv-output"] :: Bool)

	# We replace any $<parameter> pattern in `pp[save-model]` as soon as
	# possible, so if `parameter` is written wrong, the error is thrown early.
	if !isempty(pp["save-model"]) && verbose
		save_model_expanded = expand_template(pp["save-model"], pp)
		@show save_model_expanded
	end

	verbose && append!(timings, ["build_and_solve_time", "build_time"])
	instance = round_instance(instance_, pp)
	print_instance_stats(instance, verbose)

	throw_if_timeout_now(start, limit)

	solver_pp = create_unprefixed_subset(pp["solver"], pp)
	m = empty_configured_model(Val(Symbol(pp["solver"])), solver_pp)
	throw_if_timeout_now(start, limit)

	model_pp = create_unprefixed_subset(pp["model"], pp)
	bsr, mbp = build_model(
		problem, formulation, instance, m, model_pp
	)
	if verbose
		println("build_stop_reason = $(bsr)")
		println("build_stop_code = $(Int(bsr))")
	end

	if bsr == BUILT_MODEL
		if pp["relax2lp"]
			JuMP.relax_integrality(m)
		end
		num_vars = num_variables(m)
		num_constrs = num_all_constraints(m)
		if verbose
			@show num_vars
			@show num_constrs
		end
	end

	verbose && close_and_print!(timings, "build_time")
	throw_if_timeout_now(start, limit)

	if !isempty(pp["save-model"]) && verbose
		if bsr == BUILT_MODEL
			save_model_time = @elapsed begin
				@timeit TIMER "save_model" JuMP.write_to_file(m, save_model_expanded)
			end
			@show save_model_time
		else
			@warn "No model saved to '$(save_model_expanded)' because the" *
				" reason for stopping the model building was not $(BUILT_MODEL)" *
				" but instead $(bsr)."
		end
	end

	if pp["do-not-solve"]
		verbose && close_and_print!(timings, "build_and_solve_time")
		return nothing
	end

	throw_if_timeout_now(start, limit)

	if bsr == BUILT_MODEL
		verbose && println("MARK_FINAL_GENERIC_SOLVE")
		fms_name = "finished_model_solve"
		fms_time = @elapsed begin
			# The version taking an amount of seconds does not throw a timeout
			# error, it just sets the time limit.
			@timeit TIMER fms_name optimize_within_time_limit!(
				m, limit - (time() - start)
			)
		end
		if verbose
			println("$fms_name = $fms_time")
			stop_reason = termination_status(m)
			@show stop_reason
			stop_code = Int64(stop_reason)
			@show stop_code
		end
	end
	if verbose
		close_time = close_and_print!(timings, "build_and_solve_time")
		push!(timings, TimeSection("solution_print_time", close_time))
	end

	# This needs to be done even if it is not printed to warm-start the JIT
	# when a mock run is done.
	@assert bsr in (BUILT_MODEL, FOUND_OPTIMUM)
	solution = if bsr == BUILT_MODEL
		if primal_status(m) == MOI.FEASIBLE_POINT
			obj_value = objective_value(m)
			verbose && @show obj_value
			obj_bound = objective_bound(m)
			verbose && @show obj_bound
			# If the model was solved relaxed we try not to obtain a solution from
			# it (because it would not be a valid CutPattern). Also, we do not
			# print any info because we are formulation agnostic, and without
			# knowledge of the formulation, the value of the variables is of little
			# relevance.
			if pp["relax2lp"]
				nothing
			else
				get_cut_pattern(problem, formulation, m, mbp)
			end
		else
			nothing
		end
	elseif bsr == FOUND_OPTIMUM
			get_cut_pattern(problem, formulation, m, mbp)
	end
	if solution !== nothing
		@timeit TIMER "stringfy_solutions" begin
			solution_str = "solution = $solution\n"
			local pretty_sol_str :: String
			local pretty_simple_sol_str :: String
			if isa(solution, CutPattern)
				pretty_sol_str = to_pretty_str(solution)
				pretty_simple_sol_str = to_pretty_str(simplify!(deepcopy(solution)))
			else
				pretty_sol_str = join(to_pretty_str.(solution), "\n")
				pretty_simple_sol_str = join(
					to_pretty_str.(simplify!.(deepcopy.(solution))), "\n"
				)
			end

			pretty_sol_str = "PRETTY_STR_SOLUTION_BEGIN\n" * pretty_sol_str *
				"\nPRETTY_STR_SOLUTION_END\n"
			pretty_simple_sol_str = "SIMPLIFIED_PRETTY_STR_SOLUTION_BEGIN\n" *
				pretty_simple_sol_str *
				"\nSIMPLIFIED_PRETTY_STR_SOLUTION_END\n"
		end

		if verbose
			@timeit TIMER "print_solutions" begin
				iob = IOBuffer()
				write(iob, solution_str)
				write(iob, pretty_sol_str)
				write(iob, pretty_simple_sol_str)
				print(read(seekstart(iob), String))
			end
			sol_val = solution_value(problem, instance, solution)
			println("solution_value = $sol_val")
		end
	end
	verbose && close_and_print!(timings, "solution_print_time")

	# This is done just before returning because the version of
	# `optimize_within_time_limit!` called above does not throw. The
	# non-throwing version was used to be able to print as much information as
	# possible before throwing. However, if the timeout was reached, this
	# method MUST throw.
	if termination_status(m) == MOI.TIME_LIMIT
		throw(TimeoutError(start, limit, time()))
	end
	throw_if_timeout_now(start, limit)

	return
end

# Read the instance, build the model, solve the model, and print related stats.
"""
    read_build_solve_and_print(problem, format, pp)

!!! **Internal use.**

Given the parsed parameters (`pp`), read the instance file, build the model,
solve the model (unless `pp['do-not-solve']` is true), and print statistics
related to this process (unless `pp['no-csv-output']` is true).
The list of options recognized and implemented by this method is the
list returned by `generic_args()` (it also needs the required arguments
in `core_argparse_settings()`). The other arguments are solver or model
specific and are extracted and passed to their specific methods.
"""
@timeit TIMER function read_build_solve_and_print(
	format, pp
) # pp stands for parsed parameters
	verbose = !(pp["no-csv-output"] :: Bool)
	instance_path = pp["instance_path"] :: String
	verbose && @show instance_path
	timings = TimeSection[]
	verbose && push!(timings, "read_instance_time")
  dataset = Data.read_from_file(format, instance_path)
	verbose && close_and_print!(timings, "read_instance_time")
	problem = Val(Symbol(pp["problem"] :: String))
	model_type = Val(Symbol(pp["model"] :: String))

	if is_collection(format)
		for (index, instance) in enumerate(dataset)
			if verbose
				println("INSTANCE_START_MARKER_$index")
				push!(timings, "total_instance_time")
			end
			build_solve_and_print(problem, model_type, instance, pp, timings)
			if verbose
				close_and_print!(timings, "total_instance_time")
				println("INSTANCE_CLOSE_MARKER_$index")
			end
		end
	else
		if verbose
			println("INSTANCE_START_MARKER_1")
			verbose && push!(timings, "total_instance_time")
		end
		build_solve_and_print(problem, model_type, dataset, pp, timings)
		if verbose
			verbose && close_and_print!(timings, "total_instance_time")
			println("INSTANCE_CLOSE_MARKER_1")
		end
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
`solver`, and `instance_path`. They cannot be modeled as `Arg` objects because,
by design, all extra arguments must be options (i.e., be optional and preceded
by dashes).
"""
function core_argparse_settings() :: ArgParseSettings
	s = ArgParse.ArgParseSettings()
	ArgParse.add_arg_group!(s, "Core Parameters", "core_parameters")
	@add_arg_table! s begin
		"problem"
			help = "The type of problem to be solved (case sensitive, ex.: G2KP, G2CSP, G2ODP, G2MKP). Required."
			arg_type = String
		"format"
			help = "The format of the instance(s) in the file (case sensitive, ex.: Classic_G2KP, CPG_SLOPP, CPG_SSSCSP, CPG_ODPW, CPG_MHLOPPW). Required."
			arg_type = String
		"model"
			help = "Model or solution procedure to be used (case sensitive, ex.: Flow, KnapsackPlates, PPG2KP). Required."
			arg_type = String
		"solver"
			help = "Solver to be used if necessary (case-sensitive, ex.: NoSolver, Cbc, CPLEX, Gurobi). Required, even if --do-not-solve is specified. NoSolver use `JuMP.Model()` which allows building and saving the model, but not solving it."
			arg_type = String
		"instance_path"
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
			"Defines a time limit (in seconds) to be observed in the context of the model-agnostic process of reading, building, solving, and printing. Each model has to define its own flag to control the time inside the model building process (to be set independently, it is does not interact with this flag). `JuMP.set_time_limit_sec` is called over the model before starting to solve it, with the remaining time after reading instance and building the model (not the total time). If `--warm-jit` defines that there will be a mock run to warm the jit, the time limit applies to this mock run (i.e., if the mock run breaks the limit it an exception is raised) but the timer is reset after the mock, so the mock time is not counted for the 'real' run time limit. A `GuillotineModels.TimeoutError` is raised if the time limit is violated."
		),
		Arg(
			"save-model", "",
			"Save the model of the passed instance to the given filename. The format used depends on the extension of the filename, allowed formats are described by the enumeration `MathOptInterface.FileFormats.FileFormat` from package `MathOptInterface`). The filename also allows using \$<parameter_name> patterns inside the name, these will be replaced by the value of the parameter (if it was not passed, the default value will be used). For convenience, a pseudo-parameter instance_file is also available (it is the same as `basename(instance_path)`, it includes any file extensions), and `Bool` parameters (i.e., flags with no argument) are replaced by 0 or 1. Putting an @ in front of \$<some_string> will disable the substitution (and remove the @ from the final string). There is no way to express an @ followed by a substitution that actually occurs (this is a limitation of the code)."
		),
		Arg(
			"warm-jit", "no",
			"There are three valid options: 'no', 'with-toy', and 'yes'. If 'yes', the instance is solved two times, but the first time the arguments are changed to disable output (i.e., '--no-csv-output' and '--SOLVER-no-output' are passed, and if supported, '--MODEL-quiet' is passed too; in the GuillotineModels.TIMER the timings of this first run are under 'warm_jit'). If 'with-toy', it behaves similarly to 'yes' but instead of using the instance itself it uses a very small hardcoded instance (this helps a lot, but many procedures only called when the model has specific properties are not called, and therefore, not compiled). If 'no', the code is just run a single time."
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
			"round-nearest", 1.0,
			"Multiplies the instances size-related fields by the passed factor and round them to nearest (the model answer becomes a GUESS, not a valid primal heuristic, nor a valid bound)."
		),
		Arg(
			"round-up", 1.0,
			"Multiplies the instances size-related fields by the passed factor and round them up (the model becomes a PRIMAL HEURISTIC)."
		),
		Arg(
			"round-down", 1.0,
			"Multiplies the instance size-related fields by the passed factor and round them down (the model becomes an OPTIMISTIC GUESS, A VALID BOUND)."
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
	# Remember that for now, only solver and model support options/flags.
	# problem and format do not support them yet.
	core_arg_names = ["problem", "format", "solver", "model", "instance_path"]
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
		["round-nearest", "round-up", "round-down"]
	)))
	num_rounds > 1 && @error(
		"only one of --round-{nearest,up,down} may be passed at the" *
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
    toy_instance(problem, format):: String

Return a `String` representing an instance of problem in `format`.

This is the function that should be extended in order to be able to call `run`
with the `--warm-jit with-toy`.
"""
function toy_instance(problem, format)
	error(
		"If you want to mock a run with `--warm-jit with-toy` there need to" *
		" exist an implementation of `GuillotineModels.CommandLine." *
		"toy_instance` taking the instance format you passed. You have passed" *
		" problem `$problem` and format `$format`."
	)
end

function toy_instance(::Val{:G2KP}, ::Val{:CPG_SLOPP}) :: String
	return toy_instance(Val(:G2KP), CPG_SLOPP{Int,Int,Int}())
end

function toy_instance(
	::Val{:G2KP}, ::CPG_SLOPP{D, S, P}
) :: String where {D, S, P}
"""
***2D Rectangular Problem***
***Instances for the Single Large Object Placement Problem (SLOPP)***
Input parameter file: SLOPP_parameters.txt
***************************************************************************************************************
Total number of instances
LargeObject.Length      LargeObject.Width
Number of different item types (i)
Item[i].Length  Item[i].Width   Item[i].LowerBoundDemand        Item[i].UpperBoundDemand        Item[i].Value
***************************************************************************************************************
1
100	200
2
50     100    0     1     2
50     200    0     1     1
"""
end

function toy_instance(::Val{:G2MKP}, ::Val{:CPG_MHLOPPW}) :: String
	return toy_instance(Val(:G2MKP), CPG_MHLOPPW{Int, Int, Int}())
end
function toy_instance(
	::Val{:G2MKP}, ::CPG_MHLOPPW{D, S, P}
) :: String where {D, S, P}
"""
***2D Rectangular Problem***
***Instances for the Multiple Heterogeneous Large Object Placement Problem (MHLOPP/W)***
Input parameter file: written_by_hand
***************************************************************************************************************
Total number of instances
Number of different large objects (j)
LargeObject[j].Length	LargeObject[j].Width	LargeObject[j].Available
Number of different item types (i)
Item[i].Length	Item[i].Width	Item[i].LowerBoundDemand	Item[i].UpperBoundDemand	Item[i].Value
***************************************************************************************************************
1
1
100	200	2
4
100	200	0	3	20000
50	150	0	1	7501
50	200	0	2	10001
50	50	0	2	2501"""
end

const CPG_SSSCSP_HEADER = """
***2D Rectangular Problem***
***Instances for the Single Stock Size Cutting Stock Problem (SSSCSP)***
Input parameter file: written_by_hand
****************************************************************************************************
Total number of instances 
LargeObject.Length	LargeObject.Width
Number of different item types (i)
Item[i].Length	Item[i].Width	Item[i].Demand
*****************************************************************************************************"""

function toy_instance(::Val{:G2OPP}, ::Val{:CPG_SSSCSP}) :: String
	return toy_instance(Val(:G2OPP), CPG_SSSCSP{Int, Int, Int}())
end
function toy_instance(
	::Val{:G2OPP}, ::CPG_SSSCSP{D, S, P}
) :: String where {D, S, P}
"""
$CPG_SSSCSP_HEADER
1
100	200
3
100	50	2
25	100	2
50	100	1"""
end

# Same as G2OPP, but with the double of the items and therefore needing
# two original plates (i.e., large objects).
function toy_instance(::Val{:G2CSP}, ::Val{:CPG_SSSCSP}) :: String
	return toy_instance(Val(:G2CSP), CPG_SSSCSP{Int, Int, Int}())
end
function toy_instance(
	::Val{:G2CSP}, ::CPG_SSSCSP{D, S, P}
) :: String where {D, S, P}
"""
$CPG_SSSCSP_HEADER
1
100	200
3
100	50	4
25	100	4
50	100	2"""
end

function toy_instance(::Val{:G2KP}, ::Val{:Classic_G2KP}) :: String
	return toy_instance(Val(:G2KP), Classic_G2KP{Int,Int,Int}())
end

function toy_instance(
	::Val{:G2KP}, ::Classic_G2KP{D, S, P}
) :: String where {D, S, P}
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
	problem = Val(Symbol(p_args["problem"] :: String))
	format = Val(Symbol(p_args["format"] :: String))
	copyed_args = deepcopy(p_args)
	# The temporary file is only needed if p_args["warm-jit"] == "with-toy",
	# however, if it is needed, then it needs to be valid for the rest of
	# this method scope.
	mktemp() do path, io
		if p_args["warm-jit"] == "with-toy"
				write(io, toy_instance(problem, format) :: String)
				close(io)
				copyed_args["instance_path"] = path
		end
		copyed_args["no-csv-output"] = true
		copyed_args[p_args["solver"] * "-" * "no-output"] = true
		valid_model_args = accepted_arg_list(Val(Symbol(p_args["model"])))
		if any(arg -> arg.name == "quiet", valid_model_args)
			copyed_args[p_args["model"] * "-" * "quiet"] = true
		end
		read_build_solve_and_print(format, copyed_args)
	end
	return
end

# TODO: For now, there is no help support for which are possible values of
# the problem and the format positional parameters. There is also no way
# to support problem-specific (or format-specific) options/flags.
"""
    run(args = ARGS; implemented_models = [...], supported_solvers = [...])

Parse the command-line arguments and take appropriate action: prints help
or solve instance(s) in the specified file of the specified format for
the specified problem using the specified solver, model, and their options.

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
	supported_solvers :: Vector{Symbol} = [
		:CPLEX, :Gurobi, :Cbc, :GLPK, :NoSolver
	]
)
	p_args = parse_args(args, implemented_models, supported_solvers)
	isempty(p_args) && return # Happens if called just for "--help".
	verbose = !(p_args["no-csv-output"] :: Bool)
	verbose && @show args
	verbose && @show p_args
	# Out of the `if verbose` because we want it to be compiled in mock runs.
	date_now = Dates.format(Dates.now(), dateformat"yyyy-mm-ddTHH:MM:SS")
	verbose && @show date_now
	warm_jit(p_args) # It checks `p_args["warm-jit"]` to check if it will run.
	# remaining code should not depend on this option
	delete!(p_args, "warm-jit")
	# Now that is safe to extract some options, let us already call the
	# read_build_solve_and_print function with the correct format,
	# so it can be specialized for it.
	format = Val(Symbol(p_args["format"] :: String))
	total_file_time = @elapsed read_build_solve_and_print(format, p_args)
	verbose && @show total_file_time

	return
end

end # module
