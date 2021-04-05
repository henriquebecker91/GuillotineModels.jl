"""
`SaveModel` module purpose is to provide a delegate-to-solver option to
`MathOptInterface.write_to_file`. Instead of using the current MOI
fallback (that only works for MPS and is kinda memory-hungry and slow),
`SaveModel.write_to_file` checks if `SaveModel.natively_supports` returns
`true` for its arguments, and in such case calls
`SaveModel.native_write_to_file` over the object returned by the
`MOI.RawSolver()` property; otherwise it calls the `MOI` fallback.

"""
module SaveModel

import JuMP
import MathOptInterface
const MOI = MathOptInterface

import Requires.@require

import TimerOutputs.@timeit
import ..TIMER

# Based on: https://discourse.julialang.org/t/trying-to-save-mps-from-raw-solver-and-problems-with-moi-rawsolver/58515/2?u=henrique_becker
_get_raw_solver(model :: JuMP.Model) = _get_raw_solver(JuMP.backend(model))
function _get_raw_solver(model :: MOI.Utilities.CachingOptimizer)
	if model.state == MOI.Utilities.EMPTY_OPTIMIZER
		MOI.Utilities.attach_optimizer(model)
	end
	return _get_raw_solver(model.optimizer)
end
function _get_raw_solver(model :: MOI.Bridges.LazyBridgeOptimizer)
	_get_raw_solver(model.model)
end
_get_raw_solver(model) = model # Unrecognized means we arrived were we wanted.

"""
	write_to_file(model, filename; format = FORMAT_AUTOMATIC

Save the model to filename trying to use the solver mechanism instead of
the generic `JuMP`/`MathOptInterface` one.
"""
@timeit TIMER function write_to_file(
	model :: JuMP.Model,
	filename :: String;
	format :: MOI.FileFormats.FileFormat = MOI.FileFormats.FORMAT_AUTOMATIC
)
	raw_solver = _get_raw_solver(model)
	if natively_supports(raw_solver, filename; format = format)
		native_write_to_file(raw_solver, filename; format = format)
	else
		JuMP.write_to_file(model, filename; format = format)
	end

	return
end
export write_to_file

function natively_supports(
	raw_solver, filename :: String;
	format :: MOI.FileFormats.FileFormat = MOI.FileFormats.FORMAT_AUTOMATIC
) :: Bool
	return false
end

function native_write_to_file(
	raw_solver, filename :: String;
	format :: MOI.FileFormats.FileFormat = MOI.FileFormats.FORMAT_AUTOMATIC
)
	throw(ArgumentError(
		"`native_write_to_file` was called for a solver it does know." *
		" Check if the the solver package is loaded, and if some package" *
		" provides a solver-specific method for the functions in" *
		" `GuillotineModels.SaveModel`."
	))
end

const MOI_FF = MOI.FileFormats
const MOI_FF_FF = MOI_FF.FileFormat

const _KNOWN_FORMATS = Tuple{String, MOI_FF_FF}[
	(".cbf", MOI_FF.FORMAT_CBF),
	(".lp", MOI_FF.FORMAT_LP),
	(".mof.json", MOI_FF.FORMAT_MOF),
	(".mps", MOI_FF.FORMAT_MPS),
	(".sdpa", MOI_FF.FORMAT_SDPA),
	(".dat-s", MOI_FF.FORMAT_SDPA),
]

function _format_from_ext(filename :: String)
	for (ext, enum) in _KNOWN_FORMATS
		(endswith(filename, ext) || occursin("$ext.", filename)) && return enum
	end
	error("Unable to automatically detect format of $(filename).")
end

function _real_format(filename :: String, format :: MOI_FF_FF) :: MOI_FF_FF
	return if format == MOI_FF.FORMAT_AUTOMATIC
		_format_from_ext(filename) :: MOI_FF_FF
	else
		format
	end
end

function __init__()
	@require CPLEX="a076750e-1247-5638-91d2-ce28b192dca0" begin
		function natively_supports(
			model :: CPLEX.Optimizer, filename :: String;
			format :: MOI_FF_FF = MOI_FF.FORMAT_AUTOMATIC
		)
			real_format = _real_format(filename, format)
			return real_format in (MOI_FF.FORMAT_MPS, MOI_FF.FORMAT_LP)
		end

		function native_write_to_file(
			model :: CPLEX.Optimizer, filename :: String;
			format :: MOI_FF_FF = MOI_FF.FORMAT_AUTOMATIC
		)
			@assert format == MOI_FF.FORMAT_AUTOMATIC ||
				format == _format_from_ext(filename)
			return CPLEX.CPXwriteprob(model.env, model.lp, filename, C_NULL)
		end
	end # @require CPLEX

	@require Gurobi="2e9cd046-0924-5485-92f1-d5272153d98b" begin
		function natively_supports(
			model :: Gurobi.Optimizer, filename :: String;
			format :: MOI_FF_FF = MOI_FF.FORMAT_AUTOMATIC
		)
			real_format = _real_format(filename, format)
			return real_format in (MOI_FF.FORMAT_MPS, MOI_FF.FORMAT_LP)
		end
		function native_write_to_file(
			model :: Gurobi.Optimizer, filename :: String;
			format :: MOI_FF_FF = MOI_FF.FORMAT_AUTOMATIC
		)
			@assert format == MOI_FF.FORMAT_AUTOMATIC ||
				format == _format_from_ext(filename)
			# Calling GRBupdatemodel directly here is a recipe for a disaster,
			# we need to use Gurobi.jl internals to avoid any subtle bugs.
			Gurobi._update_if_necessary(model)
			return Gurobi.GRBwrite(model, filename)
		end
	end # @require Gurobi
end # __init__()

end # module
