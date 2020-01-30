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
function read_build_solve_and_print(pp) # pp stands for parsed parameters
	if !pp["no-csv-output"]
		println("instfname = $(pp["instfname"])")
	end

	N, L_, W_, l_, w_, p, d = InstanceReader.read_from_file(pp["instfname"])
	L, W, l, w = div_and_round_instance(L_, W_, l_, w_, pp)

	solver_id = Val(Symbol(pp["solver"]))
	solver_options_names = getfield.(accepted_arg_list(solver_id), :name)
	solver_pp = filter(p -> p[1] ∈ solver_options_names, pp)

	m = empty_configured_model(pp)

	model_id = Val(Symbol(pp["model"]))
	model_options_names = getfield.(accepted_arg_list(model_id), :name)
	model_pp = filter(p -> p[1] ∈ model_options_names, pp)

	# TODO: change the prefix removal to be done in a copy just before passing
	# to the relevant method (solver configuration, or model building)
	@timeit "build_model" begin
	build_model_return = build_model(
		model_id, m, d, p, l, w, L, W, model_pp
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

function gen_prefixed_argparse_settings(
	prefix :: Union{String, Symbol}
) :: ArgParse.ArgParseSettings
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

function core_argparse_settings()
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

function generic_argparse_settings()
	s = ArgParse.ArgParseSettings()
	add_arg_group(s, "Generic Options", "generic_options")
	for arg in generic_args()
		add_arg_table(s, arg)
	end
	set_default_arg_group(s)
	s
end

# NOTE: all options are parsed, even if just one model is selected at a time.
# The reasons for that are: (1) we do not know the model before the parsing
# unless we manipulate ARGS directly; (2) if they are not included in the
# parsing they do not appear in the help message.
function argparse_settings(
	models_list :: Vector{Symbol}, solvers_list :: Vector{Symbol}
)
	s = core_argparse_settings()
	ArgParse.import_settings(s, generic_argparse_settings())
	for sym in [solvers_list; models_list]
		ArgParse.import_settings(s, gen_prefixed_argparse_settings(sym))
	end
	s
end

function warn_if_unrecognized_options(p_args)
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

	valid_options_names = [
		model_arg_names; core_arg_names; generic_arg_names; solver_arg_names
	]
	for option in keys(p_args)
		if option ∉ valid_options_names
			@warn(
				"option $(option) was recognized by the parsing but is not from" *
				" the common options, nor the model selected, are you sure this" *
				" is what you wanted?"
			)
		end
	end

	return nothing
end

function throw_if_incompatible_common_options(p_args)
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
end

function remove_option_prefix(prefix, p_args)
	prefix_id = Val(Symbol(prefix))
	prefix_options_names = getfield.(accepted_arg_list(prefix_id), :name)
	for name in prefix_options_names
		prefixed_name = prefix * "-" * name
		if haskey(p_args, prefixed_name)
			p_args[name] = p_args[prefixed_name]
			delete!(p_args, prefixed_name)
		end
	end
	p_args
end

# Definition of the command line arguments.
function parse_args(args, models_list, solvers_list) :: Dict{String, Any}
	@timeit "parse_args" begin
	s = argparse_settings(models_list, solvers_list)
	p_args = ArgParse.parse_args(args, s) :: Dict{String, Any}
	warn_if_unrecognized_options(p_args)
	throw_if_incompatible_common_options(p_args)
	remove_option_prefix(p_args["model"], p_args)
	throw_if_incompatible_options(Val(Symbol(p_args["model"])), p_args)
	remove_option_prefix(p_args["solver"], p_args)
	throw_if_incompatible_options(Val(Symbol(p_args["solver"])), p_args)
	end # @timeit

	return p_args
end

function mock_instance() :: String
	"""
	100	200
	2
	50     100    1     2
	50     200    1     1
	"""
end

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
	implemented_models = [:Flow, :PPG2KP, :KnapsackPlates],
	# TODO: add the rest of the solvers below [:Cbc, :GLPK]
	supported_solvers = [:CPLEX, :Gurobi]
)
	@timeit "run" begin
	#@show args
	#@show implemented_models
	p_args = parse_args(args, implemented_models, supported_solvers)

	!p_args["do-not-mock-for-compilation"] && mock_for_compilation(p_args)
	# remaining code should not depend on this option
	delete!(p_args, "do-not-mock-for-compilation")
	total_instance_time = @elapsed read_build_solve_and_print(p_args)
	@show total_instance_time
	end # timeit
end

end # module
