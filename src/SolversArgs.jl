module SolversArgs

export accepted_arg_list, throw_if_incompatible_options
export empty_configured_model

using TimerOutputs
using JuMP
import Requires.@require

using ArgParse

import ...Utilities
import ...Utilities.Args: Arg
import ...Utilities.Args: accepted_arg_list, throw_if_incompatible_options

function Utilities.Args.accepted_arg_list(::Val{:Solvers}) :: Vector{Arg}
	return [
		Arg(
			"threads", 1,
			"Number of threads used by the solver (not model building)."
		),
		Arg(
			"time-limit", 31536000.0,
			"Time limit in seconds for solver B&B (not model building, nor root solving), the default is one year."
		),
		Arg(
			"solver-seed", 1,
			"The random seed used by the solver."
		),
		Arg(
			"no-solver-tricks", false,
			"Vary between solvers but disable things like heuristics, probe, repeat presolve, and dive."
		),
		Arg(
			"no-solver-output", false,
			"Disables the solver output."
		)
	]
end

function argparse_settings()
	s = ArgParseSettings()
	ArgParse.add_arg_group(s, "Solvers Common Flags", "solver_common_flags")
	for arg in accepted_arg_list(Val(:Solvers))
		ArgParse.add_arg_table(s, arg)
	end
	set_default_arg_group(s)
	#= TODO: if this is added back, it needs to be adapted to the new
	# style using Utilities.Args.
	ArgParse.add_arg_group(s, "CPLEX-specific Flags", "cplex_specific_flags")
	@add_arg_table s begin
		"--cplex-det-time-limit"
			help = "deterministic time limit for CPLEX B&B (not model building, nor root solving), only CPLEX has deterministic time limit among the solvers"
			arg_type = Float64
			default = [1.0E+75]
			nargs = 1
	end
	set_default_arg_group(s)
	=#
	s
end

# While solver-specific arguments and checking is not relevant, let us
# have this ::Val{Solvers} aberration.
function Utilities.Args.throw_if_incompatible_options(::Val{:Solvers}, p_args)
	# This method is called on CommandLine.jl, if conflicts arise they may
	# be just checked hered.
	#= Solver-specific flags are disabled for now
	is_default_value = p_args["cplex-det-time-limit"] == [1.0E+75]
	is_default_value && p_args["solver"] != "CPLEX" && @error(
		"flag cplex-det-time-limit was passed, but the solver selected is not CPLEX"
	)
	=#
end

function __init__()
	#function cplex_empty_configured_model(p_args)
	@require CPLEX="a076750e-1247-5638-91d2-ce28b192dca0" function empty_configured_model(::Val{:CPLEX}, p_args)
		scrind_value = p_args["no-solver-output"] ? CPLEX.CPX_OFF : CPLEX.CPX_ON
		# TODO: just append options if no-solver-tricks is on, do not
		# write the common options two times
		#=
		if p_args["no-solver-tricks"]
			configuration = [
				"CPX_PARAM_EPGAP" => 1e-6
				, "CPX_PARAM_PROBE" => -1
				, "CPX_PARAM_HEURFREQ" => -1
				, "CPX_PARAM_REPEATPRESOLVE" => -1
				, "CPX_PARAM_DIVETYPE" => 1
				#, "CPX_PARAM_DETTILIM" => p_args["cplex-det-time-limit"],
				, "CPX_PARAM_TILIM" => p_args["time-limit"]
				#, "CPX_PARAM_VARSEL" => CPLEX.CPX_VARSEL_MAXINFEAS,
				, "CPX_PARAM_STARTALG" => CPLEX.CPX_ALG_BARRIER
				, "CPX_PARAM_SCRIND" => scrind_value
				, "CPX_PARAM_THREADS" => p_args["threads"]
				, "CPX_PARAM_RANDOMSEED" => p_args["solver-seed"]
			]
		else
			configuration = [
				"CPX_PARAM_EPGAP" => 1e-6,
				#, "CPLEX.CPX_PARAM_DETTILIM" = p_args["cplex-det-time-limit"],
				"CPX_PARAM_TILIM" => p_args["time-limit"],
				"CPX_PARAM_VARSEL" => CPLEX.CPX_VARSEL_MAXINFEAS,
				"CPX_PARAM_STARTALG" => CPLEX.CPX_ALG_BARRIER,
				"CPX_PARAM_SCRIND" => scrind_value,
				"CPX_PARAM_THREADS" => p_args["threads"],
				"CPX_PARAM_RANDOMSEED" => p_args["solver-seed"]
			]
		end=#
		#=
		CPLEX.Optimizer(
			CPX_PARAM_EPGAP = 1e-6
			, CPX_PARAM_PROBE = -1
			, CPX_PARAM_HEURFREQ = -1
			, CPX_PARAM_REPEATPRESOLVE = -1
			, CPX_PARAM_DIVETYPE = 1
			#, CPX_PARAM_DETTILIM = p_args["cplex-det-time-limit"],
			, CPX_PARAM_TILIM = p_args["time-limit"]
			#, CPX_PARAM_VARSEL = CPLEX.CPX_VARSEL_MAXINFEAS,
			, CPX_PARAM_STARTALG = CPLEX.CPX_ALG_BARRIER
			, CPX_PARAM_SCRIND = scrind_value
			, CPX_PARAM_THREADS = p_args["threads"]
			, CPX_PARAM_RANDOMSEED = p_args["solver-seed"]
		)
		=#

		configuration = [
			"CPX_PARAM_EPGAP" => 1e-6,
			#, "CPLEX.CPX_PARAM_DETTILIM" = p_args["cplex-det-time-limit"],
			"CPX_PARAM_TILIM" => p_args["time-limit"],
			"CPX_PARAM_VARSEL" => CPLEX.CPX_VARSEL_MAXINFEAS,
			"CPX_PARAM_STARTALG" => CPLEX.CPX_ALG_BARRIER,
			"CPX_PARAM_SCRIND" => scrind_value,
			"CPX_PARAM_THREADS" => p_args["threads"],
			"CPX_PARAM_RANDOMSEED" => p_args["solver-seed"]
		]
		model = JuMP.direct_model(CPLEX.Optimizer())
		for c in configuration
			JuMP.set_parameter(model, c...)
		end
		model
	end

	@require Gurobi="2e9cd046-0924-5485-92f1-d5272153d98b" function empty_configured_model(::Val{:Gurobi}, p_args)
		JuMP.direct_model(
			if p_args["no-solver-tricks"]
				Gurobi.Optimizer(
					Method = 2, # use barrier for LP
					#PreSparsify = 1, # try to reduce nonzeros
					#Presolve = 2, # aggressive presolving
					Threads = p_args["threads"],
					Seed = p_args["solver-seed"],
					OutputFlag = no_solver_out ? 0 : 1,
					MIPGap = 1e-6,
					TimeLimit = p_args["time-limit"]
				)
			else
				Gurobi.Optimizer(
					Method = 2, # use barrier for LP
					#PreSparsify = 1, # try to reduce nonzeros
					#Presolve = 2, # aggressive presolving
					Threads = p_args["threads"],
					Seed = p_args["solver-seed"],
					OutputFlag = no_solver_out ? 0 : 1,
					MIPGap = 1e-6,
					TimeLimit = p_args["time-limit"]
				)
			end
		)
	end

	@require Cbc="9961bab8-2fa3-5c5a-9d89-47fab24efd76" function empty_configured_model(::Val{:Cbc}, p_args)
		JuMP.direct_model(
			# TODO: the options of cbc were not studied so, for now it does
			# the same thing independent of the value of disable-solver-tricks
			Cbc.Optimizer(
				threads = p_args["threads"],
				ratioGap = 1e-6,
				logLevel = p_args["no-solver-output"] ? 0 : 1,
				randomSeed = p_args["solver-seed"],
				barrier = true,
				seconds = p_args["time-limit"]
			)
		)
	end
end

function empty_configured_model(
	::Val{T}, p_args
) where {T}
	@error(
		"solver " * string(T) * " is not implemented, define an" *
		" implementation of empty_configured_model for it"
	)
end

# Create a new model with a configured underlying solver.
function empty_configured_model(p_args)
	@timeit "empty_configured_model" begin
	model = empty_configured_model(Val(Symbol(p_args["solver"])), p_args)
	end # timeit

	return model
end

end # module
