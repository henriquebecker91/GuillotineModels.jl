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

	function Utilities.Args.accepted_arg_list(::Val{:CPLEX}) :: Vector{Arg}
		return [
			Arg(
				"threads", 1,
				"Number of threads used by the solver (not model building)."
			),
			Arg(
				"time-limit", 31536000.0,
				"Time limit in seconds for solver B&B (not root solving), the default is one year."
			),
			Arg(
				"solver-seed", 1,
				"The random seed used by the solver."
			),
			Arg(
				"no-tricks", false,
				"TODO: rewrite this help to say specifically what is disabled."
			),
			Arg(
				"no-output", false,
				"Disables the solver output."
			)
		]
	end

	function Utilities.Args.throw_if_incompatible_options(::Val{:CPLEX}, p_args)
		# TODO: should we throw if the option does not exist?
	end
	end # @require CPLEX

	@require Gurobi="2e9cd046-0924-5485-92f1-d5272153d98b" function empty_configured_model(::Val{:Gurobi}, p_args)
		JuMP.direct_model(
			if p_args["no-solver-tricks"]
				Gurobi.Optimizer(
					Method = 2, # use barrier for LP
					#PreSparsify = 1, # try to reduce nonzeros
					#Presolve = 2, # aggressive presolving
					Threads = p_args["threads"],
					Seed = p_args["solver-seed"],
					OutputFlag = p_args["no-output"] ? 0 : 1,
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
					OutputFlag = p_args["no-output"] ? 0 : 1,
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
				logLevel = p_args["no-output"] ? 0 : 1,
				randomSeed = p_args["solver-seed"],
				barrier = true,
				seconds = p_args["time-limit"]
			)
		)
	end

	@require GLPK="60bf3e95-4087-53dc-ae20-288a0d20c6a6" function empty_configured_model(::Val{:GLPK}, p_args)
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
