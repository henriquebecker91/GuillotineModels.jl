module Args

import ...Utilities
import ...Utilities.Args: Arg
export accepted_arg_list, throw_if_incompatible_options

function Utilities.Args.accepted_arg_list(::Val{:PPG2KP})
	return [
		Arg(
			"lower-bound", 0.0,
			"The use of these values is dependent on the other options selected, check code."
		),
		Arg(
			"upper-bound", Inf,
			"The use of these values is dependent on the other options selected, check code."
		),
		Arg(
			"round2disc", false,
			"Round the second child size to a discretized position."
		),
		# TODO: edit the below after the pricing is implemented.
		Arg(
			"faithful2furini2016", false,
			"Tries to be the most faithful possible to the description on the Furini2016 paper (more specifically the complete model, with the reductions, FOR NOW without a heuristic solution first and without pricing); the flags --no-cut-position, --no-redundant-cut and --no-furini-symmbreak disable parts of this reimplementation."
		),
		Arg(
			"no-redundant-cut", false,
			"Disables Furini2016 Redundant-Cut reduction: a bunch of flags is used to check if some trim cut is necessary for optimality, or is dominated by other trim cuts."
		),
		Arg(
			"no-cut-position", false,
			"Disables Furini2016 Cut-Position reduction: if a trivial heuristic shows that no combination of six or more pieces fit a plate, then the plate may be cut with restricted cuts without loss to optimality."
		),
		Arg(
			"no-furini-symmbreak", false,
			"Disables Furini2016 symmetry breaking: as trim cuts are needed in faithful2furini2016 mode, the symmetry-breaking is less restrictive than 'do not cut after midplate', it just removes each cut y > midplate in which plate_size - y is an existing cut; enabling this flag makes the code just use all discretized positions as cuts."
		),
		Arg(
			"no-pricing", false,
			"Disables the pricing procedure (both the iterative and the final). Combined with faithful2furini2016 this will give the Complete PP-G2KP Model with Cut-Position and Redundant-Cut reductions (pass no-redundant-cut and no-cut-position to have the pure Complete PP-G2KP Model)."
		),
		#=Arg(
			"warm-start", false,
			"(For now conflict with --faithful2furini2016) Uses the heuristic described in Furini2016 to generate a initial primal feasible solution, warm-start the model with unrestricted cuts fixed to zero, and then unfix the unrestricted cuts to solve the complete model."
		),=#
		Arg(
			"ignore-d", false,
			"Ignore the demand information during discretization, used to measure impact."
		),
		Arg(
			"use-c25", false,
			"Add the tightening constraints 2.5. Increase considerably the number of constraints hoping to improve the relaxation a little."
		),
		Arg(
			"ignore-2th-dim", false,
			"Ignore the dimension not being discretized during discretization, used to measure impact."
		),
		Arg(
			"pricing-heuristic-seed", 0,
			"Defines the seed used to start the RNG object passed to the `GuillotineModels.PPG2KP.iterated_greedy` method, which result is used in the final pricing of the restricted model which is, finally, used for the final pricing of the complete model. If you pass the no-pricing flag, this is not used."
		)
	]
end

function Utilities.Args.throw_if_incompatible_options(::Val{:PPG2KP}, p_args)
end

end # module

