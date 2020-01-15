module SolverArgParse

export empty_configured_model

# TODO: remove after putting require to work
import Cbc

function get_arg_parse_settings()
	s = ArgParseSettings()
	@add_arg_table s begin
		"--solver"
			help = "which Solver to use if the model is solved (case-insensitive)"
			arg_type = String
			default = ["Cbc"]
			nargs = 1
		"--threads"
			help = "number of threads used by the solver (not model building)"
			arg_type = Int
			default = [1]
			nargs = 1
		"--time-limit"
			help = "time limit in seconds for solver B&B (not model building, nor root solving), the default is one year"
			arg_type = Float64
			default = [31536000.0]
			nargs = 1
		"--cplex-det-time-limit"
			help = "deterministic time limit for CPLEX B&B (not model building, nor root solving), only CPLEX has deterministic time limit among the solvers"
			arg_type = Float64
			default = [1.0E+75]
			nargs = 1
		"--solver-seed"
			help = "the random seed used by the solver (the model building is deterministic)"
			arg_type = Int
			default = [1]
			nargs = 1
		"--disable-solver-tricks"
			help = "vary between solvers but disable things like heuristics, probe, repeat presolve, and dive"
			nargs = 0
		"--disable-solver-output"
			help = "disable the solver output"
			nargs = 0
	end
	s
end

# TODO: the variables of these methods are not declared (they come from
# the hash). Change them back to hash accesses.
function empty_configured_model(
	::Val{:cbc}, p_args; no_solver_out = no_solver_out
)
	Base.require(:Cbc)
	JuMP.direct_model(
		# TODO: the options of cbc were not studied so, for now it does
		# the same thing independent of the value of disable-solver-tricks
		if p_args["disable-solver-tricks"]
			Cbc.Optimizer(
				threads = 1,
				ratioGap = 1e-6,
				logLevel = no_solver_out ? 0 : 1,
				randomSeed = first(p_args["solver-seed"]),
				barrier = true,
				seconds = first(p_args["time-limit"])
			)
		else
			Cbc.Optimizer(
				threads = 1,
				ratioGap = 1e-6,
				logLevel = no_solver_out ? 0 : 1,
				randomSeed = first(p_args["solver-seed"]),
				barrier = true,
				seconds = first(p_args["time-limit"])
			)
		end
	)
end

function empty_configured_model(
	::Val{:cplex}, p_args; no_solver_out = no_solver_out
)
	Base.require(:CPLEX)
	JuMP.direct_model(
		if p_args["disable-solver-tricks"]
			CPLEX.Optimizer(
				CPX_PARAM_EPGAP = 1e-6,
				CPX_PARAM_PROBE = -1,
				CPX_PARAM_HEURFREQ = -1,
				CPX_PARAM_REPEATPRESOLVE = -1,
				CPX_PARAM_DIVETYPE = 1,
				CPX_PARAM_DETTILIM = first(p_args["det-time-limit"]),
				CPX_PARAM_TILIM = first(p_args["time-limit"]),
				#CPX_PARAM_VARSEL = CPLEX.CPX_VARSEL_MAXINFEAS,
				CPX_PARAM_STARTALG = CPLEX.CPX_ALG_BARRIER,
				CPX_PARAM_SCRIND = no_solver_out ? CPLEX.CPX_OFF : CPLEX.CPX_ON,
				CPX_PARAM_THREADS = first(p_args["threads"]),
				CPX_PARAM_RANDOMSEED = first(p_args["solver-seed"])
			)
		else
			CPLEX.Optimizer(
				CPX_PARAM_EPGAP = 1e-6,
				CPX_PARAM_DETTILIM = first(p_args["det-time-limit"]),
				CPX_PARAM_TILIM = first(p_args["time-limit"]),
				CPX_PARAM_VARSEL = CPLEX.CPX_VARSEL_MAXINFEAS,
				CPX_PARAM_STARTALG = CPLEX.CPX_ALG_BARRIER,
				CPX_PARAM_SCRIND = no_solver_out ? CPLEX.CPX_OFF : CPLEX.CPX_ON,
				CPX_PARAM_THREADS = first(p_args["threads"]),
				CPX_PARAM_RANDOMSEED = first(p_args["solver-seed"])
			)
		end
	)
end

function empty_configured_model(
	::Val{:gurobi}, p_args; no_solver_out = no_solver_out
)
	Base.require(:Gurobi)
	JuMP.direct_model(
		if p_args["disable-solver-tricks"]
			Gurobi.Optimizer(
				Method = 2, # use barrier for LP
				#PreSparsify = 1, # try to reduce nonzeros
				#Presolve = 2, # aggressive presolving
				Threads = first(p_args["threads"]),
				Seed = first(p_args["solver-seed"]),
				OutputFlag = no_solver_out ? 0 : 1,
				MIPGap = 1e-6,
				TimeLimit = first(p_args["time-limit"])
			)
		else
			Gurobi.Optimizer(
				Method = 2, # use barrier for LP
				#PreSparsify = 1, # try to reduce nonzeros
				#Presolve = 2, # aggressive presolving
				Threads = first(p_args["threads"]),
				Seed = first(p_args["solver-seed"]),
				OutputFlag = no_solver_out ? 0 : 1,
				MIPGap = 1e-6,
				TimeLimit = first(p_args["time-limit"])
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
	solver_id = Val(Symbol(lowercase(first(p_args["solver"]))))
	return empty_configured_model(
		solver_id, p_args; no_solver_out = no_solver_out
	)
end

end # module
