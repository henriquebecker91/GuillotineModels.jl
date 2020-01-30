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
		configuration = [
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
			)
		]
	end

	function Utilities.Args.throw_if_incompatible_options(::Val{:CPLEX}, p_args)
		# TODO: should we throw if the option does not exist?
	end
	end # @require CPLEX

	@require Gurobi="2e9cd046-0924-5485-92f1-d5272153d98b" begin
	function empty_configured_model(::Val{:Gurobi}, p_args)
		# Note: without the typeassert below some values are casted and Gurobi.jl
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
			)
		]
	end

	function Utilities.Args.throw_if_incompatible_options(::Val{:Gurobi}, p_args)
		# TODO: should we throw if the option does not exist?
	end
	end # @require Gurobi

	@require Cbc="9961bab8-2fa3-5c5a-9d89-47fab24efd76" begin
	function empty_configured_model(::Val{:Cbc}, p_args)
		JuMP.direct_model(
			# TODO: the options of cbc were not studied so, for now it does
			# the same thing independent of the value of disable-solver-tricks
			Cbc.Optimizer(
				threads = p_args["threads"],
				ratioGap = 1e-6,
				logLevel = p_args["no-output"] ? 0 : 1,
				randomSeed = p_args["seed"],
				barrier = true,
				seconds = p_args["time-limit"]
			)
		)
	end
	end # @require Cbc

	@require GLPK="60bf3e95-4087-53dc-ae20-288a0d20c6a6" begin
	function empty_configured_model(::Val{:GLPK}, p_args)
		model = JuMP.direct_model(GLPK.Optimizer())
		base_conf = [
			# GLPK tm_lim is in milisseconds
			"tm_lim" => p_args["time-limit"] * 1000,
			"msg_lev" => (p_args["no-output"] ? GLPK.MSG_OFF : GLPK.MSG_ON)
		]
		for c in configuration
			JuMP.set_parameter(model, c...)
		end
		model
	end
	end # @require GLPK
end

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

# Create a new model with a configured underlying solver.
function empty_configured_model(p_args)
	@timeit "empty_configured_model" begin
	model = empty_configured_model(Val(Symbol(p_args["solver"])), p_args)
	end # timeit

	return model
end

end # module
