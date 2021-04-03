"""
`SolversArgs` module aggregates the methods necessary for making solvers
available to use with `GuillotineModels.CommandLine` methods (especially
`run`).

The package `Requires.jl` is used to implement
`GuillotineModels.SolversArgs.empty_configured_model(Val(:SOME_SOLVER))`
only if `SOME_SOLVER` package is loaded.

The solver-specific methods this module implements for each supported solver
are: `GuillotineModels.Utilities.Args.{accepted_arg_list,
throw_if_incompatible_options}`, and `GuillotineModels.CommandLine.SolversArgs.
empty_configured_model`. The `GuillotineModels.Utilities` methods are always
available, not just when the respective solver package was imported.

The supported solvers are `CPLEX`, `Gurobi`, `Cbc`, and `GLPK`.
Supporting new solvers is just a question of implementing a specialization
of the three methods mentioned above for the specific solver (see the
code for examples).
"""
module SolversArgs

export empty_configured_model

import TimerOutputs.@timeit
import ...TIMER
import Requires.@require
import JuMP

using ArgParse

import ...Utilities
import ...Utilities.Args: Arg
import ...Utilities.Args: accepted_arg_list, throw_if_incompatible_options

# CPLEX is not define in this scope, so we return a symbol with the
# name of the constant and use `getfield` inside the scope CPLEX
# is defined.
const CPLEX_NAME2CODE_LP_ALG = Dict{String, Symbol}(
	"CPX_ALG_AUTOMATIC" => :CPX_ALG_AUTOMATIC,
	"auto" => :CPX_ALG_AUTOMATIC,
	"automatic" => :CPX_ALG_AUTOMATIC,
	"0" => :CPX_ALG_AUTOMATIC,
	"CPX_ALG_PRIMAL" => :CPX_ALG_PRIMAL,
	"primal" => :CPX_ALG_PRIMAL,
	"1" => :CPX_ALG_PRIMAL,
	"CPX_ALG_DUAL" => :CPX_ALG_DUAL,
	"dual" => :CPX_ALG_DUAL,
	"2" => :CPX_ALG_DUAL,
	"CPX_ALG_NET" => :CPX_ALG_NET,
	"net" => :CPX_ALG_NET,
	"network" => :CPX_ALG_NET,
	"3" => :CPX_ALG_NET,
	"CPX_ALG_BARRIER" => :CPX_ALG_BARRIER,
	"barrier" => :CPX_ALG_BARRIER,
	"4" => :CPX_ALG_BARRIER,
	"CPX_ALG_SIFTING" => :CPX_ALG_SIFTING,
	"sifting" => :CPX_ALG_SIFTING,
	"5" => :CPX_ALG_SIFTING,
	"CPX_ALG_CONCURRENT" => :CPX_ALG_CONCURRENT,
	"concurrent" => :CPX_ALG_CONCURRENT,
	"6" => :CPX_ALG_CONCURRENT,
	"CPX_ALG_CONCURRENT" => :CPX_ALG_CONCURRENT,
	"parallel" => :CPX_ALG_CONCURRENT,
	"7" => :CPX_ALG_CONCURRENT
)
function __init__()
	#function cplex_empty_configured_model(p_args)
	@require CPLEX="a076750e-1247-5638-91d2-ce28b192dca0" begin
	@timeit TIMER function empty_configured_model(::Val{:CPLEX}, p_args)
		scrind_value = p_args["no-output"] ? CPLEX.CPX_OFF : CPLEX.CPX_ON
		root_relax_method = CPLEX_NAME2CODE_LP_ALG[p_args["root-relax-method"]]
		lp_method = CPLEX_NAME2CODE_LP_ALG[p_args["LP-method"]]
		# https://www.ibm.com/support/pages/cplex-performance-tuning-linear-programs
		configuration = Pair{String, Any}[
			"CPX_PARAM_EPGAP" => 1e-6,
			"CPX_PARAM_TILIM" => p_args["time-limit"],
			"CPX_PARAM_STARTALG" => getfield(CPLEX, root_relax_method),
			"CPX_PARAM_LPMETHOD" => getfield(CPLEX, lp_method),
			# "Sifting is a simple form of column generation well suited for models
			# where the number of variables dramatically exceeds the number of
			# constraints."
			#"CPX_PARAM_STARTALG" => CPLEX.CPX_ALG_SIFTING,
			#"CPX_PARAM_LPMETHOD" => CPLEX.CPX_ALG_SIFTING,
			# "The barrier method tends to work well on problems where the product
			# of the constraint matrix multiplied by its transpose is sparse. "
			#"CPX_PARAM_STARTALG" => CPLEX.CPX_ALG_BARRIER,
			#"CPX_PARAM_LPMETHOD" => CPLEX.CPX_ALG_BARRIER,
			#"CPX_PARAM_STARTALG" => CPLEX.CPX_ALG_NET,
			#"CPX_PARAM_LPMETHOD" => CPLEX.CPX_ALG_NET,
			#"CPX_PARAM_BARDISPLAY" => 2, # 2 == diagnostic information level
			#"CPX_PARAM_SIMDISPLAY" => 2, # 2 == diagnostic information level
			# For the LPs of the iterative pricing of PPG2KP we need to avoid
			# numerical instability problems.
			#"CPX_PARAM_BARALG" => 1,
			# "[...] the computation time for the simplex method depends more on
			# the number of constraints than the number of variables."
			# And the dual is the opposite (CPLEX.CPX_ALG_DUAL).
			#"CPX_PARAM_STARTALG" => CPLEX.CPX_ALG_PRIMAL,
			#"CPX_PARAM_LPMETHOD" => CPLEX.CPX_ALG_PRIMAL,
			#"CPX_PARAM_PPRIIND" => CPLEX.CPX_PPRIIND_FULL, # pricing inside simplex
			#"CPX_PARAM_PERIND" => CPLEX.CPX_ON, # start using perturbations
			#"CPX_PARAM_PERLIM" => 1000000, # num degenerate iters until perturbation
			# Group parameter to help with numerical instability without the need
			# of fine-tuning.
			#"CPX_PARAM_NUMERICALEMPHASIS" => CPLEX.CPX_ON,
			#"CPX_PARAM_BAREPCOMP" => 1e-10,
			#"CPX_PARAM_EPMRK" => 0.9, # last measure against numerical instability
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
		JuMP.set_optimizer_attributes(model, configuration...)
		model
	end
	end # @require CPLEX

	@require Gurobi="2e9cd046-0924-5485-92f1-d5272153d98b" begin
	@timeit TIMER function empty_configured_model(::Val{:Gurobi}, p_args)
		# Note: without the type prefix below some values are casted and Gurobi.jl
		# fails to find the correct method to call (it differs from Int and Float).
		configuration = Pair{String, Any}[
			"Method" => p_args["LP-method"],
			# "NodeMethod" => 2, # use barrier for MIP non-root nodes
			"Threads" => p_args["threads"],
			"OutputFlag" => p_args["no-output"] ? 0 : 1,
			"TimeLimit" => p_args["time-limit"],
			"Seed" => p_args["seed"],
			"MIPGap" => 1e-6
		]
		# TODO: implement a no-trick flag for Gurobi?
		raw_parameters = eval(
			Meta.parse(p_args["raw-parameters"])
		) :: Vector{Pair{String, Any}}
		configuration = [configuration; raw_parameters]
		model = JuMP.direct_model(Gurobi.Optimizer(Gurobi.Env()))
		JuMP.set_optimizer_attributes(model, configuration...)
		model
	end
	end # @require Gurobi

	@require Cbc="9961bab8-2fa3-5c5a-9d89-47fab24efd76" begin
	@timeit TIMER function empty_configured_model(::Val{:Cbc}, p_args)
		configuration = Pair{String, Any}[
			# TODO: discover if the barrier algorithm is being correctly used.
			# The `barrier` command-line "parameter" is in fact an action/command
			# that solves the model relaxation with the barrier algorithm, and
			# not an argument for which linear algorithm the branch-and-bound
			# should use in the root node. The `logLevel` parameter, for example,
			# does not work with barrier, which will output even if `logLevel` is
			# set to zero.
			"logLevel" => p_args["no-output"] ? 0 : 1
			, "threads" => p_args["threads"]
			, "randomCbcSeed" => p_args["seed"]
			, "seconds" => p_args["time-limit"]
			, "ratioGap" => 1e-6
			, "barrier" => ""
		]
		# TODO: implement a no-trick flag for Cbc?
		raw_parameters = eval(
			Meta.parse(p_args["raw-parameters"])
		) :: Vector{Pair{String, Any}}
		configuration = [configuration; raw_parameters]
		model = JuMP.Model(JuMP.with_optimizer(Cbc.Optimizer))
		JuMP.set_optimizer_attributes(model, configuration...)
		model
	end
	end # @require Cbc

	@require GLPK="60bf3e95-4087-53dc-ae20-288a0d20c6a6" begin
	@timeit TIMER function empty_configured_model(::Val{:GLPK}, p_args)
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
		JuMP.set_optimizer_attributes(model, configuration...)
		model
	end
	end # @require GLPK
end

function Utilities.Args.accepted_arg_list(::Val{:CPLEX}) :: Vector{Arg}
	return [
		Arg(
			"root-relax-method", "automatic",
			"Changes CPX_PARAM_LPMETHOD. Accept simplified single-word names (primal, dual, etc...), CPLEX names (CPX_ALG_PRIMAL, CPX_ALG_DUAL, ...), and codes (1, 2, ...)."
		),
		Arg(
			"LP-method", "automatic",
			"Changes CPX_PARAM_STARTALG/CPXPARAM_MIP_Strategy_StartAlgorithm. Accept simplified single-word names (primal, dual, etc...), CPLEX names (CPX_ALG_PRIMAL, CPX_ALG_DUAL, ...), and codes (1, 2, ...)."
		),
		Arg(
			"threads", 1,
			"The value of CPXPARAM_Threads. If zero and no callbacks is the number os cores, if zero and callbacks is one. If a positive number, is that number of cores."
		),
		Arg(
			"time-limit", 31536000.0,
			"BROKEN, DO NOT USE, ALWAYS OVERWRITTEN BY `--generic-time-limit`. The value of CPX_PARAM_TILIM. Depends on CPX_PARAM_CLOCKTYPE which is by default wall clock time (not CPU time). Time limit in seconds for solver B&B (not root solving). Our default is one year."
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
	warn_about_time_limit(Val(:CPLEX), p_args)
	valid_lp_alg_names = sort!(collect(keys(CPLEX_NAME2CODE_LP_ALG)))
	lpm = "LP-method"
	rrm = "root-relax-method"
	Utilities.throw_if_unrecognized(lpm, p_args[lpm], valid_lp_alg_names)
	Utilities.throw_if_unrecognized(rrm, p_args[rrm], valid_lp_alg_names)
end

function Utilities.Args.accepted_arg_list(::Val{:Gurobi}) :: Vector{Arg}
	return [
		Arg(
			"LP-method", -1,
			"Gurobi has a parameter called 'Method' that defines the algorithm used to solve continuous models (including MIP root node continuous relaxations). The options (described in Gurobi 9.0 reference) are: automatic (-1), primal simplex (0), dual simplex (1), barrier (2), concurrent (3), deterministic concurrent (4), deterministic concurrent simplex (5)."
		),
		Arg(
			"threads", 1,
			"Number of threads for all Gurobi parallelizable algorithms. Zero is automatic, probably the number of cores but may be fewer. If a positive number, is that number of cores."
		),
		Arg(
			"time-limit", 31536000.0,
			"BROKEN, DO NOT USE, ALWAYS OVERWRITTEN BY `--generic-time-limit`. Set Gurobi parameter: TimeLimit. Total time limit in seconds. Our default is one year."
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
	warn_about_time_limit(Val(:Gurobi), p_args)
	# TODO: should we throw if the option does not exist?
end

function Utilities.Args.accepted_arg_list(::Val{:Cbc}) :: Vector{Arg}
	return [
		Arg(
			"threads", 1,
			"Number of threads for \"parallel branch-and-bound\"."
		),
		Arg(
			"time-limit", 31536000.0,
			"BROKEN, DO NOT USE, ALWAYS OVERWRITTEN BY `--generic-time-limit`. Set Cbc parameter: seconds. Total time limit in seconds (? not very well documented)."
		),
		Arg(
			"seed", 1,
			"Set Cbc parameter: randomCbcSeed. Our default is 1 (different from Cbc default that is -1)."
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
	warn_about_time_limit(Val(:Cbc), p_args)
	# TODO: should we throw if the option does not exist?
end

function Utilities.Args.accepted_arg_list(::Val{:GLPK}) :: Vector{Arg}
	return [
		Arg(
			"time-limit", 2097152,
			"BROKEN, DO NOT USE, ALWAYS OVERWRITTEN BY `--generic-time-limit`. Set GLPK parameter: tm_lim. The original parameter is in milliseconds, but to keep it similar to the other solvers this option takes seconds. To set this parameter with milliseconds precision use the --raw-parameter option. The default is a little over 24 days because GLPK uses a Int32 for milliseconds and 24 days is close to the maximum time-limit supported."
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

function warn_about_time_limit(::Val{SOLVER_SYM}, p_args) where {SOLVER_SYM}
	if haskey(p_args, "time-limit")
		tl = p_args["time-limit"]
		args = Utilities.Args.accepted_arg_list(Val(SOLVER_SYM))
		tl_arg_idx = findfirst(arg -> arg.name == "time-limit", args)
		if !isnothing(tl_arg_idx) && args[tl_arg_idx].default != tl
			@warn "--$SOLVER_SYM-time-limit is BROKEN, it is always overwritten" *
				" by --generic-time-limit even when this one is not passed (i.e., it" *
				" is overwritten by its default value in such case)."
		end
	end

	return
end

function Utilities.Args.throw_if_incompatible_options(::Val{:GLPK}, p_args)
	warn_about_time_limit(Val(:GLPK), p_args)
	# TODO: should we throw if the option does not exist?
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
		"Solver " * string(T) * " is not supported (i.e., there is not" *
		" an implementation of Julia method `GuillotineModels.CommandLine" *
		".SolversArgs.empty_configured_model(::Val{Symbol(\"$(string(T))\")}," *
		"p_args)` for it) or the solver package was not imported before this" *
		"method was called."
	)
end

function empty_configured_model(
	::Val{:NoSolver}, p_args
) where {T}
	return JuMP.Model()
end

function Utilities.Args.accepted_arg_list(::Val{:NoSolver}) :: Vector{Arg}
	return Arg[
		Arg(
			"no-output", false,
			"Provided to avoid breaking the assumption every solver provides this" *
			" flag. Does notthing, as NoSolver should not print anything anyway."
		),
	]
end

function Utilities.Args.throw_if_incompatible_options(::Val{:NoSolver}, p_args)
	warn_about_time_limit(Val(:NoSolver), p_args)
	return nothing
end

end # module
