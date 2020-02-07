"""
SolversArgs module aggregates the methods necessary for making solvers
available to use with GuillotineModels.CommandLine methods.

The package Requires is used to implement solver-specific methods only if the
solver package is loaded. The only non-solver-specific method implemented
is the error fallback for `empty_configured_model`.

The solver-specific methods this module implements for each supported solver
are: `GuillotineModels.Utilities.Args.{accepted_arg_list,
throw_if_incompatible_options}`, and `GuillotineModels.CommandLine.SolversArgs.
empty_configured_model`.

The supported solvers are `CPLEX`, `Gurobi`, `Cbc`, and `GLPK`.
Supporting new solvers is just a question of implementing a specialization
of the three methods mentioned above for the specific solver (see the
code for examples).
"""
module SolversArgs

export empty_configured_model

import TimerOutputs # needed for the expansion of @timeit
import TimerOutputs.@timeit
import Requires.@require
import JuMP

using ArgParse

import ...Utilities
import ...Utilities.Args: Arg
import ...Utilities.Args: accepted_arg_list, throw_if_incompatible_options

function __init__()
	#function cplex_empty_configured_model(p_args)
	@require CPLEX="a076750e-1247-5638-91d2-ce28b192dca0" begin
	function empty_configured_model(::Val{:CPLEX}, p_args)
		scrind_value = p_args["no-output"] ? CPLEX.CPX_OFF : CPLEX.CPX_ON
		configuration = Pair{String, Any}[
			"CPX_PARAM_EPGAP" => 1e-6,
			#, "CPLEX.CPX_PARAM_DETTILIM" = p_args["cplex-det-time-limit"],
			"CPX_PARAM_TILIM" => p_args["time-limit"],
			"CPX_PARAM_STARTALG" => CPLEX.CPX_ALG_BARRIER,
			"CPX_PARAM_SCRIND" => scrind_value,
			"CPX_PARAM_THREADS" => p_args["threads"],
			"CPX_PARAM_RANDOMSEED" => p_args["seed"]
		]
		if p_args["no-tricks"]
			extra_conf = [
				"CPX_PARAM_PROBE" => -1
				, "CPX_PARAM_HEURFREQ" => -1
				, "CPX_PARAM_REPEATPRESOLVE" => -1
				, "CPX_PARAM_DIVETYPE" => 1
			]
			configuration = [configuration; extra_conf]
		end
		raw_parameters = eval(
			Meta.parse(p_args["raw-parameters"])
		) :: Vector{Pair{String, Any}}
		configuration = [configuration; raw_parameters]
		model = JuMP.direct_model(CPLEX.Optimizer())
		for c in configuration
			JuMP.set_parameter(model, c...)
		end
		model
	end

	function Utilities.Args.accepted_arg_list(::Val{:CPLEX}) :: Vector{Arg}
		return [
			Arg(
				"threads", 1,
				"The value of CPXPARAM_Threads. If zero and no callbacks is the number os cores, if zero and callbacks is one. If a positive number, is that number of cores."
			),
			Arg(
				"time-limit", 31536000.0,
				"The value of CPX_PARAM_TILIM. Depends on CPX_PARAM_CLOCKTYPE which is by default wall clock time (not CPU time). Time limit in seconds for solver B&B (not root solving). Our default is one year."
			),
			Arg(
				"seed", 1,
				"The value of CPX_PARAM_RANDOMSEED. The random seed used by CPLEX."
			),
			Arg(
				"no-tricks", false,
				"Set CPX_PARAM_PROBE, CPX_PARAM_HEURFREQ, and CPX_PARAM_REPEATPRESOLVE to -1. Also, set CPX_PARAM_DIVETYPE to 1. Basically, disable many tricks used by CPLEX to speedup the solving (but that sometimes have the opposite effect)."
			),
			Arg(
				"no-output", false,
				"Set CPX_PARAM_SCRIND to false. Disables the solver output."
			),
			Arg(
				"raw-parameters", "Pair{String, Any}[]",
				"A string of Julia code to be evaluated to CPLEX parameters. For example: 'Pair{String, Any}[\"CPX_PARAM_SCRIND\" => CPLEX.CPX_OFF]' will have the same effect as --CPLEX-no-output. Obviously, this is a security breach."
			)
		]
	end

	function Utilities.Args.throw_if_incompatible_options(::Val{:CPLEX}, p_args)
		# TODO: should we throw if the option does not exist?
	end
	end # @require CPLEX

	@require Gurobi="2e9cd046-0924-5485-92f1-d5272153d98b" begin
	function empty_configured_model(::Val{:Gurobi}, p_args)
		# Note: without the type prefix below some values are casted and Gurobi.jl
		# fails to find the correct method to call (it differs from Int and Float).
		configuration = Pair{String, Any}[
			"Method" => 2 # use barrier for LP
			, "Threads" => p_args["threads"]
			, "OutputFlag" => p_args["no-output"] ? 0 : 1
			, "TimeLimit" => p_args["time-limit"]
			, "Seed" => p_args["seed"]
			, "MIPGap" => 1e-6
		]
		# TODO: implement a no-trick flag for Gurobi?
		raw_parameters = eval(
			Meta.parse(p_args["raw-parameters"])
		) :: Vector{Pair{String, Any}}
		configuration = [configuration; raw_parameters]
		model = JuMP.direct_model(Gurobi.Optimizer())
		for c in configuration
			JuMP.set_parameter(model, c...)
		end
		model
	end

	function Utilities.Args.accepted_arg_list(::Val{:Gurobi}) :: Vector{Arg}
		return [
			Arg(
				"threads", 1,
				"Number of threads for all Gurobi parallelizable algorithms. Zero is automatic, probably the number of cores but may be fewer. If a positive number, is that number of cores."
			),
			Arg(
				"time-limit", 31536000.0,
				"Set Gurobi parameter: TimeLimit. Total time limit in seconds. Our default is one year."
			),
			Arg(
				"seed", 1,
				"The random seed used by Gurobi. Our default (one) is different from Gurobi (which is zero)."
			),
			Arg(
				"no-output", false,
				"Set OutputFlag to zero. Disables the solver output."
			),
			Arg(
				"raw-parameters", "Pair{String, Any}[]",
				"A string of Julia code to be evaluated to Gurobi parameters. For example: 'Pair{String, Any}[\"OutputFlag\" => 0]' will have the same effect as --Gurobi-no-output. Obviously, this is a security breach."
			)

		]
	end

	function Utilities.Args.throw_if_incompatible_options(::Val{:Gurobi}, p_args)
		# TODO: should we throw if the option does not exist?
	end
	end # @require Gurobi

	@require Cbc="9961bab8-2fa3-5c5a-9d89-47fab24efd76" begin
	function empty_configured_model(::Val{:Cbc}, p_args)
		configuration = Pair{String, Any}[
			"barrier" => true # use barrier for LP
			, "threads" => p_args["threads"]
			, "logLevel" => p_args["no-output"] ? 0 : 1
			, "seconds" => p_args["time-limit"]
			, "randomSeed" => p_args["seed"]
			, "ratioGap" => 1e-6
		]
		# TODO: implement a no-trick flag for Cbc?
		raw_parameters = eval(
			Meta.parse(p_args["raw-parameters"])
		) :: Vector{Pair{String, Any}}
		configuration = [configuration; raw_parameters]
		model = JuMP.Model(JuMP.with_optimizer(Cbc.Optimizer))
		for c in configuration
			JuMP.set_parameter(model, c...)
		end
		model
	end

	function Utilities.Args.accepted_arg_list(::Val{:Cbc}) :: Vector{Arg}
		return [
			Arg(
				"threads", 1,
				"Number of threads for \"parallel branch-and-bound\"."
			),
			Arg(
				"time-limit", 31536000.0,
				"Set Cbc parameter: seconds. Total time limit in seconds (? not very well documented)."
			),
			Arg(
				"seed", 1,
				"Set Cbc parameter: randomSeed. Our default is 1 (different from Cbc default that is 1234567)."
			),
			Arg(
				"no-output", false,
				"Set logLevel to zero. Disables the solver output."
			),
			Arg(
				"raw-parameters", "Pair{String, Any}[]",
				"A string of Julia code to be evaluated to Cbc parameters. For example: 'Pair{String, Any}[\"logLevel\" => 0]' will have the same effect as --Cbc-no-output. Obviously, this is a security breach. The complete list of parameters can be found by running the cbc executable and typing ? at the prompt."
			)

		]
	end

	function Utilities.Args.throw_if_incompatible_options(::Val{:Cbc}, p_args)
		# TODO: should we throw if the option does not exist?
	end
	end # @require Cbc

	@require GLPK="60bf3e95-4087-53dc-ae20-288a0d20c6a6" begin
	function empty_configured_model(::Val{:GLPK}, p_args)
		configuration = Pair{String, Any}[
			# TODO: find how configure: # of threads, tolerance gap, barrier, etc...
			# TODO: specially: does GLPK have a RANDOMSEED option?!
			"msg_lev" => p_args["no-output"] ? GLPK.MSG_OFF : GLPK.MSG_ON
			, "tm_lim" => p_args["time-limit"] * 1000 # is milliseconds not seconds
		]
		# TODO: implement a no-trick flag for GLPK? does GLPK even have tricks?
		raw_parameters = eval(
			Meta.parse(p_args["raw-parameters"])
		) :: Vector{Pair{String, Any}}
		configuration = [configuration; raw_parameters]
			model = JuMP.direct_model(GLPK.Optimizer())
		for c in configuration
			JuMP.set_parameter(model, c...)
		end
		model
	end

	function Utilities.Args.accepted_arg_list(::Val{:GLPK}) :: Vector{Arg}
		return [
			Arg(
				"time-limit", 2097152,
				"Set GLPK parameter: tm_lim. The original parameter is in milliseconds, but to keep it similar to the other solvers this option takes seconds. To set this parameter with milliseconds precision use the --raw-parameter option. The default is a little over 24 days because GLPK uses a Int32 for milliseconds and 24 days is close to the maximum time-limit supported."
			),
			Arg(
				"no-output", false,
				"Set msg_lev to GLPK.MSG_OFF. Disables the solver output."
			),
			Arg(
				"raw-parameters", "Pair{String, Any}[]",
				"A string of Julia code to be evaluated to GLPK parameters. For example: 'Pair{String, Any}[\"msg_lev\" => GLPK.MSG_OFF]' will have the same effect as --GLPK-no-output. Obviously, this is a security breach. If you find the complete list of GLPK parameters please send it to me (henriquebecker91@gmail.com)."
			)

		]
	end

	function Utilities.Args.throw_if_incompatible_options(::Val{:GLPK}, p_args)
		# TODO: should we throw if the option does not exist?
	end
	end # @require GLPK
end

"""
    empty_configured_model(::Val{T}, p_args) :: JuMP-like Model

Creates an empty (no variables or constraints) and configured (`p_args` is used
to set the parameters of an attached solver) model of the solver `::Val{T}`
(for example, passing `::Val{:CPLEX}` will call the specialized method for
CPLEX solver).

The options available for `p_args` are available in `GuillotineModels.
Utilities.Args.accepted_arg_list(::Val{SOLVER_PACKAGE_SYMBOL})`.
"""
function empty_configured_model(
	::Val{T}, p_args
) where {T}
	@error(
		"solver " * string(T) * " is not supported (i.e., there is not" *
		" an implementation of Julia method 'empty_configured_model(" *
		"::Val{:SolverPackageName}, p_args)' for it) or the solver package" *
		" was not imported before this method was called."
	)
end

end # module
