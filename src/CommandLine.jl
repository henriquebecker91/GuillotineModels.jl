module CommandLine

# External packages used.
using TimerOutputs
using JuMP
using Random # for rng object passed to heuristic
import ArgParse
import ArgParse: set_default_arg_group, add_arg_group, @add_arg_table
# Both MathOptInterface and MathOptFormat are used to allow saving the
# generated model.
using MathOptInterface
const MOI = MathOptInterface
using MathOptFormat

include("./SolversArgs.jl") # empty_configured_model, *parse_settings()
using .SolversArgs

using ..InstanceReader # for reading the instances
using ..Utilities # for useful helper methods not implemented in JuMP

# The implemented models.
import ..Flow
import ..PPG2KP
import ..KnapsackPlates

# Style guideline: as the module block is left unindented, the @timeit
# blocks that wrap the whole method body also are not indented.

function save_model(model, filename = "saved_model.mps") :: Nothing
	@timeit "save_model" begin
	mps_model = MathOptFormat.MPS.Model()
	MOI.copy_to(mps_model, backend(model))
	MOI.write_to_file(mps_model)
	end # timeit
	nothing
end

function div_and_round_instance(L, W, l, w, p_args)
	@timeit "div_and_round_instance" begin
	# assert explanation: at least two of the three flags are disabled (i.e., 
	# have value one)
	@assert sum(isone.((
		p_args["div-and-round-nearest"][1],
		p_args["div-and-round-up"][1],
		p_args["div-and-round-down"][1]
	))) >= 2

	if p_args["div-and-round-nearest"][1] != 1
		factor = p_args["div-and-round-nearest"][1]
		roundmode = RoundNearest
	elseif p_args["div-and-round-up"][1] != 1
		factor = p_args["div-and-round-up"][1]
		roundmode = RoundUp
	elseif p_args["div-and-round-down"][1] != 1
		factor = p_args["div-and-round-down"][1]
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

function get_submodule(module_name)
	@show module_name
	@show @__MODULE__
	@show parentmodule(@__MODULE__)
	@show names(parentmodule(@__MODULE__), all = true)
	f = getfield(parentmodule(@__MODULE__), Symbol(module_name))
	@show f
	f
end

function get_submodule_field(module_name, method_name)
	submodule = get_submodule(module_name)
	getfield(submodule, Symbol(method_name))
end

# Read the instance, build the model, solve the model, and print related stats.
function read_build_solve_and_print(
	p_args, instfname_idx; no_csv_out = false, no_solver_out = false
)
	instfname = p_args["instfnames"][instfname_idx]

	if !no_csv_out
		@show instfname
		seed = first(p_args["solver-seed"])
		@show seed
	end

	L_, W_, l_, w_, p, d = InstanceReader.read_instance(instfname)
	L, W, l, w = div_and_round_instance(L_, W_, l_, w_, p_args)

	m = empty_configured_model(p_args; no_solver_out = no_solver_out)

	@timeit "build_model" begin
	model_builder = get_submodule_field(p_args["model"], :build)
	model_build_output = model_builder(
		m, d, p, l, w, L, W; p_args = p_args
	)
	end
	time_to_build_model = TimerOutputs.time(get_defaulttimer(), "build_model")

	p_args["save-model"] && !no_csv_out &&
		save_model(m, "./$(basename(instfname)).mps")

	if !no_csv_out
		@show time_to_build_model
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
		if !p_args["flow-model"] && !p_args["break-hvcut-symmetry"]
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

	if p_args["warm-start"]
		@assert !p_args["faithful2furini2016"] && !p_args["flow-model"]
		(heur_obj, heur_sel, heur_pat), heur_time, _, _, _ = @timed iterated_greedy(
			d, p, l, w, L, W, MersenneTwister(p_args["solver-seed"][1])
		)
		if !no_csv_out
			@show heur_time
			@show heur_obj
			@show heur_sel
			@show heur_pat
		end
		saved_bounds = PPG2KP.disable_unrestricted_cuts!(
			m, sort(l), sort(w), hvcuts, pli_lwb
		)
		heur_ws_time = @elapsed PPG2KP.warm_start(
			m, l, w, L, W, heur_pat, pli_lwb, hvcuts, np;
			faithful2furini2016 = p_args["faithful2furini2016"]
		)
		!no_csv_out && @show heur_ws_time
		restricted_time = @elapsed optimize!(m)
		!no_csv_out && @show restricted_time
		restricted_sol = value.(all_variables(m))
		restricted_objval = objective_value(m)
		PPG2KP.restore_bound!.(saved_bounds)
		set_start_value.(all_variables(m), restricted_sol)
	end

	p_args["do-not-solve"] && return nothing

	vars_before_deletes = all_variables(m)
	if p_args["final-pricing"] || p_args["relax2lp"]
		original_settings = relax_all_vars!(m)
	end
	time_to_solve_model = @elapsed optimize!(m)
	if (p_args["final-pricing"] || p_args["relax2lp"]) && !no_csv_out
		println("time_to_solve_relaxed_model = $(time_to_solve_model)")
	end
	if p_args["final-pricing"]
		# get the given lower bound on the instance if there is one
		best_lb = get(p_args["lower-bounds"], instfname_idx, 0)
		# if warm-start is enabled and the lb is better, use it
		p_args["warm-start"] && (best_lb = max(best_lb, restricted_objval))
		# delete all variables which would only reduce obj to below the lower bound
		#if !no_csv_out
		#  @profile (mask = delete_vars_by_pricing!(m, Float64(best_lb)))
		#  Profile.print()
		#else
		#  mask = delete_vars_by_pricing!(m, Float64(best_lb))
		#end
		kept = which_vars_to_keep(m, Float64(best_lb))
		unrelax_vars!(vars_before_deletes, m, original_settings)
		saved_bounds = fix_vars!(all_variables(m)[.!kept])
		if !no_csv_out
			num_vars_after_pricing = sum(kept)
			num_vars_del_by_pricing = length(kept) - num_vars_after_pricing
			@show num_vars_after_pricing
			@show num_vars_del_by_pricing
		end
		time_to_solve_model += @elapsed optimize!(m)
	end
	sleep(0.01) # shamefully needed to avoid out of order messages from cplex

	if !no_csv_out
		@show time_to_solve_model
		if primal_status(m) == MOI.FEASIBLE_POINT
			obj_value = objective_value(m)
		else
			obj_value = 0
		end
		@show obj_value
		if !p_args["relax2lp"]
			obj_bound = objective_bound(m)
			@show obj_bound
		end
		stop_reason = termination_status(m)
		@show stop_reason
		stop_code = Int64(stop_reason)
		@show stop_code
		if p_args["flow-model"]
			ps = value.(m[:edge][pii] for pii = 1:length(d) if is_valid(m, m[:edge][pii]))
			ps_nz = [iv for iv in enumerate(ps) if iv[2] > 0.001]
			@show ps_nz
		else
			if p_args["break-hvcut-symmetry"]
				ps = m[:pieces_sold]
				cm = m[:cuts_made]
				ps_nz = [(i, value(ps[i])) for i = 1:length(ps) if is_valid(m, ps[i]) && value(ps[i]) > 0.001]
				cm_nz = [(i, value(cm[i])) for i = 1:length(cm) if is_valid(m, cm[i]) && value(cm[i]) > 0.001]
				@show ps_nz
				@show cm_nz
				println("(piece length, piece width) => (piece index, amount in solution, profit of single piece, total profit contributed in solution) ")
				foreach(ps_nz) do e
					pii, v = e
					println((l[pii], w[pii]) => (pii, v, p[pii], v * p[pii]))
				end
				println("(parent plate length, parent plate width) => ((first child plate length, first child plate width), (second child plate length, second child plate width))")
				foreach(cm_nz) do e
					i, _ = e
					parent, fchild, schild = hvcuts[i]
					if iszero(schild)
						println((pli2lwsb[parent][1], pli2lwsb[parent][2]) => ((pli2lwsb[fchild][1], pli2lwsb[fchild][2]), (0, 0)))
					else
						println((pli2lwsb[parent][1], pli2lwsb[parent][2]) => ((pli2lwsb[fchild][1], pli2lwsb[fchild][2]), (pli2lwsb[schild][1], pli2lwsb[schild][2])))
					end
				end
			else # same as: if !p_args["break-hvcut-symmetry"]
				p_args["relax2lp"] && @warn "relax2lp flag used, be careful regarding solution"
				ps_nz = Vector{Tuple{Int, Float64}}()
				cm_nz = Vector{Tuple{Int, Float64}}()
				for i = 1:length(ps)
					if is_valid(m, ps[i]) && value(ps[i]) > 0.0001
						push!(ps_nz, (i, value(ps[i])))
					end
				end
				for i = 1:length(cm)
					if is_valid(m, cm[i]) && value(cm[i]) > 0.0001
						push!(cm_nz, (i, value(cm[i])))
					end
				end
				#ps_nz = [(i, value(ps[i])) for i = 1:length(ps) if (is_valid(m, ps[i]) && value(ps[i]) > 0.001)]
				#cm_nz = [(i, value(cm[i])) for i = 1:length(cm) if (is_valid(m, cm[i]) && value(cm[i]) > 0.001)]
				@show ps_nz
				@show cm_nz
				println("(plate length, plate width) => (number of times this extraction happened, piece length, piece width)")
				foreach(ps_nz) do e
					i, v = e
					pli, pii = np[i]
					println((pli_lwb[pli][1], pli_lwb[pli][2]) => (v, l[pii], w[pii]))
				end
				println("(parent plate length, parent plate width) => (number of times this cut happened, (first child plate length, first child plate width), (second child plate length, second child plate width))")
				foreach(cm_nz) do e
					i, v = e
					parent, fchild, schild = hvcuts[i]
		@assert !iszero(parent)
		@assert !iszero(fchild)
		if iszero(schild)
			@assert p_args["faithful2furini2016"]
			if pli_lwb[parent][1] == pli_lwb[fchild][1]
							println((pli_lwb[parent][1], pli_lwb[parent][2]) => (v, (pli_lwb[fchild][1], pli_lwb[fchild][2]), (pli_lwb[parent][1], pli_lwb[parent][2] - pli_lwb[fchild][2])))
			else
				@assert pli_lwb[parent][2] == pli_lwb[fchild][2]
							println((pli_lwb[parent][1], pli_lwb[parent][2]) => (v, (pli_lwb[fchild][1], pli_lwb[fchild][2]), (pli_lwb[parent][1] - pli_lwb[fchild][1], pli_lwb[parent][2])))
			end
		else
						println((pli_lwb[parent][1], pli_lwb[parent][2]) => (v, (pli_lwb[fchild][1], pli_lwb[fchild][2]), (pli_lwb[schild][1], pli_lwb[schild][2])))
		end
				end
			end
		end
	end

	return nothing
end

function core_parse_settings()
	s = ArgParse.ArgParseSettings()
	ArgParse.add_arg_group(s, "Core Parameters", "core_parameters")
	@add_arg_table s begin
		"model"
			help = "Model or solution procedure to be used (case sensitive, ex.: Flow, KnapsackPlates, PPG2KP). Required."
			arg_type = String
		"solver"
			help = "Solver to be used if necessary (case-sensitive, ex.: Cbc, CPLEX, Gurobi). Required, even if --do-not-solve is specified."
			arg_type = String
		"instfnames"
			help = "The paths to the instances to be solved."
			#action = :store_arg
			nargs = '*' # can pass no instances to call "-h"
	end
	set_default_arg_group(s)
	s
end

function generic_parse_settings()
	s = ArgParse.ArgParseSettings()
	add_arg_group(s, "Generic Flags", "generic_flags")
	@add_arg_table s begin
		"--do-not-solve"
			help = "The model is build but not solved. A solver has yet to be specified. Note that just building a model may depend on using a solver over subproblems. Such uses of the solver are not disabled by this flag."
			nargs = 0
		"--save-model"
			help = "Save the model of each problem instance to the working directory ('./instance_name.mps'). Uses MPS format."
			nargs = 0
		"--do-not-mock-first-for-compilation"
			help = "To avoid counting the compilation time, the first instance is solved twice, this flag disables the first and silent solve."
			nargs = 0
		"--relax2lp"
			help = "Integer and binary variables become continuous."
			nargs = 0
		"--lower-bounds"
			help = "Takes a single string, the string has N comma-separated-numbers where N is the number of instances passed, no spaces allowed, it is up to the model generation procedure to use this option."
			nargs = 1
		"--upper-bounds"
			help = "Takes a single string, the string has N comma-separated-numbers where N is the number of instances passed, no spaces allowed, it is up to the model generation procedure to use this option."
			nargs = 1
		"--div-and-round-nearest"
			help = "Divide the instances lenght and width (but not profit) by the passed factor and round them to nearest (the model answer becomes a GUESS, not a valid primal heuristic, nor a valid bound)."
			arg_type = Int
			default = [1]
			nargs = 1
		"--div-and-round-up"
			help = "Divide the instances lenght and width (but not profit) by the passed factor and round them up (the model becomes a PRIMAL HEURISTIC)."
			arg_type = Int
			default = [1]
			nargs = 1
		"--div-and-round-down"
			help = "Divide the instances lenght and width (but not profit) by the passed factor and round them down (the model becomes an OPTIMISTIC GUESS, A VALID BOUND)."
			arg_type = Int
			default = [1]
			nargs = 1
	end
	set_default_arg_group(s)
	s
end

function common_parse_settings()
	s = core_parse_settings()
	ArgParse.import_settings(s, generic_parse_settings())
	ArgParse.import_settings(s, SolversArgs.common_parse_settings())
	s
end

# NOTE: all flags are parsed, even if just one model is selected at a time.
# The reasons for that are: (1) we do not know the model before the parsing
# unless we manipulate ARGS directly; (2) if they are not included in the
# parsing they do not appear in the help message.
function parse_settings()
	s = common_parse_settings()
	ArgParse.import_settings(s, PPG2KP.Args.parse_settings())
	ArgParse.import_settings(s, Flow.Args.parse_settings())
	ArgParse.import_settings(s, KnapsackPlates.Args.parse_settings())
	s
end

function model_parse_settings(
	model_module_name :: Union{AbstractString,Symbol}
)
	args_submod = get_submodule_field(model_module_name, :Args)
	s = getfield(args_submod, :parse_settings)()
end

function extract_option_names(s)
	# NOTE: these fields are not described in the documentation, so this is
	# probably not guaranteed against sudden internal changes in the ArgParse
	# package.
	getfield.(s.args_table.fields, :dest_name)
end

function check_flag_conflicts(p_args)
	# Generic Flags conflicts
	num_rounds = sum(.! isone.(map(flag -> p_args[flag][1],
		["div-and-round-nearest", "div-and-round-up", "div-and-round-down"]
	)))
	num_rounds > 1 && @error(
		"only one of --div-and-round-{nearest,up,down} may be passed at the" *
		" same time (what the fuck you expected?)"
	)
	# Check if all options are from common or the model selected.
	# NOTE: if someday we have solver-specific arguments we need to treat them
	# here too.
	valid_option_names = [extract_option_names(common_parse_settings());
		extract_option_names(model_parse_settings(p_args["model"]))]
	for option in keys(p_args)
		if option ∉ valid_option_names
			@warn(
				"option $(option) was recognized by the parsing but is not from" *
				" the common options, nor the model selected, are you sure this" *
				" is what you wanted?"
			)
		end
	end
	# Model specific conflicts
	p_args["model"] == "PPG2KP" && PPG2KP.Args.check_flag_conflicts(p_args)
	p_args["model"] == "Flow" && Flow.Args.check_flag_conflicts(p_args)
	p_args["model"] == "KnapsackPlates" &&
		KnapsackPlates.Args.check_flag_conflicts(p_args)
end

# TODO: this could be parsed with custom type parsing:
# https://carlobaldassi.github.io/ArgParse.jl/stable/custom.html
# However, the size checks would need to go to other place.
function bounds_extra_parse(p_args)
	for option in ["lower-bounds", "upper-bounds"]
		if !isempty(p_args[option])
			@assert isone(length(p_args[option]))
			str_list = p_args[option][1]
			if match(r"^([1-9][0-9]*|0)(,([1-9][0-9]*|0))*$", str_list) === nothing
				@error("the --$(option) does not follow the desired format:" *
					" integer,integer,integer,... (no comma after last number)"
				)
			end
			num_list = parse.(Int64, split(str_list, ','))
			if length(num_list) != length(p_args["instfnames"])
				@error("the length of the list passed to --$(option) should be" *
					" the exact same size as the number of instances provided"
				)
			end
			p_args[option] = num_list
		end
	end
end

# Definition of the command line arguments.
function parse_args(args = ARGS)
	@timeit "parse_args" begin
	s = parse_settings()
	p_args = ArgParse.parse_args(args, s)
	bounds_extra_parse(p_args)
	check_flag_conflicts(p_args)
	end # timeit

	return p_args
end

# Parse the command line arguments, and call the solve for each instance.
function run(args = ARGS)
	@timeit "run" begin
	@show args
	p_args = parse_args(args)

	mock_first = !p_args["do-not-mock-first-for-compilation"]
	if mock_first && !isempty(p_args["instfnames"])
		@timeit "mock" begin
		read_build_solve_and_print(
			p_args, 1; no_csv_out = true, no_solver_out = true
		)
		end
	end
	for inst_idx in eachindex(p_args["instfnames"])
		total_instance_time = @elapsed read_build_solve_and_print(
			p_args, inst_idx; no_solver_out = p_args["disable-solver-output"]
		)
		@show total_instance_time
	end
	end # timeit
end

end # module
