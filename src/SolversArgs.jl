module SolversArgs

export empty_configured_model

using TimerOutputs
using ArgParse
using JuMP
#import CPLEX

function common_parse_settings()
	s = ArgParseSettings()
	ArgParse.add_arg_group(s, "Solvers Common Flags", "solver_common_flags")
	@add_arg_table s begin
		"--threads"
			help = "Number of threads used by the solver (not model building)."
			arg_type = Int
			default = 1
			action = :store_arg
		"--time-limit"
			help = "Time limit in seconds for solver B&B (not model building, nor root solving), the default is one year."
			arg_type = Float64
			default = 31536000.0
			action = :store_arg
		"--solver-seed"
			help = "The random seed used by the solver."
			arg_type = Int
			default = 1
			action = :store_arg
		"--disable-solver-tricks"
			help = "Vary between solvers but disable things like heuristics, probe, repeat presolve, and dive."
			action = :store_true
		"--disable-solver-output"
			help = "Disables the solver output."
			action = :store_true
	end
	set_default_arg_group(s)
	s
end

function parse_settings()
	#=
	ArgParse.add_arg_group(s, "CPLEX-specific Flags", "cplex_specific_flags")
	@add_arg_table s begin
		"--cplex-det-time-limit"
			help = "deterministic time limit for CPLEX B&B (not model building, nor root solving), only CPLEX has deterministic time limit among the solvers"
			arg_type = Float64
			default = [1.0E+75]
			nargs = 1
	end
	=#
	set_default_arg_group(s)
	s
end

function check_flag_conflicts(p_args)
	# Some error checking.
	#= Solver-specific flags are disabled for now
	is_default_value = p_args["cplex-det-time-limit"] == [1.0E+75]
	is_default_value && p_args["solver"] != "CPLEX" && @error(
		"flag cplex-det-time-limit was passed, but the solver selected is not CPLEX"
	)
	=#
end

function empty_configured_model(
	::Val{:Cbc}, p_args; no_solver_out = no_solver_out
)
	Base.require(:Cbc)
	JuMP.direct_model(
		# TODO: the options of cbc were not studied so, for now it does
		# the same thing independent of the value of disable-solver-tricks
		Cbc.Optimizer(
			threads = p_args["threads"],
			ratioGap = 1e-6,
			logLevel = no_solver_out ? 0 : 1,
			randomSeed = p_args["solver-seed"],
			barrier = true,
			seconds = p_args["time-limit"]
		)
	)
end

function cplex_inner_call(p_args, no_solver_out)
	scrind_value = no_solver_out ? CPLEX.CPX_OFF : CPLEX.CPX_ON
	# TODO: just append options if disable-solver-tricks is on, do not
	# write the common options two times
	#=
	if p_args["disable-solver-tricks"]
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

function empty_configured_model(
	::Val{:CPLEX}, p_args; no_solver_out = no_solver_out
)
	#Base.require(
	#	Base.PkgId(Base.UUID("a076750e-1247-5638-91d2-ce28b192dca0"), "CPLEX")
	#)
	Base.include_string(@__MODULE__, "import CPLEX")
	Base.invokelatest(cplex_inner_call, p_args, no_solver_out)
	#CPLEX = getfield(Main, :CPLEX)
	#:CPLEX)
end

function empty_configured_model(
	::Val{:Gurobi}, p_args; no_solver_out = no_solver_out
)
	Base.require(:Gurobi)
	JuMP.direct_model(
		if p_args["disable-solver-tricks"]
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

function empty_configured_model(
	::Val{T}, p_args; no_solver_out = no_solver_out
) where {T}
	@error("model " * T * "is not implemented, define an implementation of" *
		" empty_configured_model for it")
end

# Create a new model with a configured underlying solver.
function empty_configured_model(p_args; no_solver_out = no_solver_out)
	@timeit "empty_configured_model" begin
	solver_sym = Symbol(p_args["solver"])

	@show solver_sym
	model = empty_configured_model(
		Val(solver_sym), p_args; no_solver_out = no_solver_out
	)
	@show model
	end # timeit

	return model
end

end # module
