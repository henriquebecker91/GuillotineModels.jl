module Args

import ...Utilities
import ...Utilities.Args: Arg
export accepted_arg_list, throw_if_incompatible_options

function Utilities.Args.accepted_arg_list(::Val{:PPG2KP})
	return [
		Arg(
			"minimal-subset-in-iterative-pricing", false,
			"Only used if the value of the 'pricing' option is 'furini'. The explanation of the iterative pricing algorithm in 10.1287/ijoc.2016.0710 implies all rows are present since start, and that cuts (`x` variables) that cannot assume non-zero values yet (because the plates they are cutting are not made available by other cuts yet) are also present since start. This flag enables an alternative behavior, in which the only plate constraints in the initial model are either: obtainable by restricted cuts (i.e., present in a restricted version of the model); or have a sell/extraction (`y` variable in original model or `e` variable in the revised model) selling that plate (or extracting a piece from it). All selling/extraction variables are present, but the cut variables are a smaller subset (only the cuts over 'reachable' plates, no cuts over plates that cannot be yet generated by the model). Consequently, both 'cut variables'/columns and 'plate constraints'/rows are added dinamically/iteratively to the model."
		),
		Arg(
			"no-sort-in-iterative-pricing", false,
			"Only used if the value of the 'pricing' option is 'furini'. The explanation of the pricing algorithm in 10.1287/ijoc.2016.0710 does not make very clear if the variables with positive reduced profit added each iteration are the ones with highest reduced profit or an arbitrary subset. The default is to sort because: the extra effort is minimal compared to solving the LP; the time solving the LP can be drastically reduced by the adding the variables with highest reduced profit first."
		),
		Arg(
			"building-time-limit", 60.0 * 60.0 * 24.0 * 365.0,
			"Defines a time limit (in seconds) to be observed in the context of the model building proccess. If the solver is called during the building procedure then `JuMP.set_time_limit_sec` is called over the model before starting to solve it, with the remaining time (not the total time). A `GuillotineModels.TimeoutError` is raised if the time limit is violated."
		),
		Arg(
			"Gurobi-LP-method-inside-furini-pricing", -2,
			"Allows to switch the algorithm used to solve continuous models and root node relaxations inside the Furini's iterative pricing when Gurobi is the solver used. Will have no effect but print warnings if either Furini's pricing is not called or Gurobi is not the solver in use. The default value of -2 will not touch Gurobi's Method parameter, values -1 to 5 will set Gurobi's Method to the corresponding value during the subprocedure mentioned, and restore the previous value after it."
		),
		Arg(
			"pricing-alpha", 0.20,
			"Used to compute the number of variables added in each iteration of the iterative pricing (only used if furini pricing is selected, see pricing flag). Must be above zero and at most one. Explained in 10.1287/ijoc.2016.0710, p. 13 (747) (last paragraph before section 4.3). Default value is the one used in the experiments of the original paper."
		),
		Arg(
			"pricing-beta", 0.25,
			"Used to compute the number of variables addded in each iteration of the iterative pricing (only used if furini pricing is selected, see pricing flag). Must be non-negative, and makes most sense to be at most one. Explained in 10.1287/ijoc.2016.0710, p. 13 (747) (last paragraph before section 4.3). Default value is the one used in the experiments of the original paper."
		),
		Arg(
			"round2disc", false,
			"Round the second child size to a discretized position."
		),
		Arg(
			"hybridize-with-restricted", false,
			"If a cut matches the length or width of a piece, and this length/width cannot be obtained by any combination of smaller pieces, then the cut already immediatelly extracts the matching piece. Should preserve optimality but it is experimental."
		),
		Arg(
			"faithful2furini2016", false,
			"Tries to be the most faithful possible to the description on the Furini2016 paper of the PPG2KP (i.e., the model with the CutPosition and RedundantCut reductions BUT NOT the multistep pricing procedure); the flags --no-cut-position, --no-redundant-cut, --no-furini-symmbreak, can disable parts of this reimplementation. Passing '--pricing furini' will enable the pricing making the model the Priced PPG2KP."
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
			"pricing", "none",
			"Define if a pricing procedure will be used, and which one. The default value is 'none'. The accepted values are: 'furini', 'becker', and 'none'. Combining '--PPG2KP-pricing none' with '--PPG2KP-faithful2furini2016' will give the Complete PP-G2KP Model with Cut-Position and Redundant-Cut reductions (pass no-redundant-cut and no-cut-position to have the pure Complete PP-G2KP Model)."
		),
		Arg(
			"ignore-d", false,
			"Ignore the demand information during discretization, used to measure impact. The flag `hybridize-with-restricted` is affected by it too (i.e., less hybridizations will happen)."
		),
		Arg(
			"use-c25", false,
			"Add the tightening constraints 2.5. Increase considerably the number of constraints hoping to improve the relaxation a little."
		),
		Arg(
			"ignore-2th-dim", false,
			"Ignore the dimension not being discretized during discretization, used to measure impact. The flag `hybridize-with-restricted` is affected by it too (i.e., less hybridizations will happen)."
		),
		Arg(
			"heuristic-seed", 1,
			"Used only if the value of the 'MIP-start' option implies the use of `GuillotineModels.PPG2KP.Heuristic.fast_iterated_greedy` (the seed is used to initialize a RNG, more specifically a `RandomNumbers.Xoroshiro128Plus`)."
		),
		Arg(
			"verbose", false,
			"Only really relevant if you are studying the code. Will increase" *
			" the amount of data outputted. The quiet flag overrides this."
		),
		Arg(
			"quiet", false,
			"Avoid outputting any information or warning, just output errors." *
			" Overrides the verbose flag. ."
		),
		Arg(
			"do-not-purge-unreachable", false,
			"By default, the code removes variables and constraints that are" *
			" useless (because the only way to reach them was removed by some" *
			" other mechanism, like pricing), this flag disables this feature" *
			" and leave the work to the presolver of the MIP solver. Note that" *
			" if the 'verbose' flag is enabled then the unreachable variables" *
			" and constraints will be searched and printed, even if they are" *
			" not removed from the model after."
		),
		Arg("MIP-start", "expected",
			"Three values are allowed: \"expected\", \"none\", and \"guaranteed\". If \"expected\" is passed, then the model will not be MIP-started unless some sort of pricing (that already generates a complete feasible solution in the process) is used. If it is Becker's pricing, then the model is just MIP-started with the simple iterated greedy 2-staged heuristic; if it is Furini's pricing, then two MIP-starts occur: the first one using the same heuristic in a restricted model; the second one using the solution of the restricted model in the returned non-restricted model. If \"none\" is passed, then the model is never MIP-started, even temporarily (as the first MIP-start that can occur in the middle of Furini's pricing). If \"guaranteed\" is passed, it behaves as \"expected\" if pricing is enabled, but if pricing is disabled, it calls the iterated greedy 2-staged heuristic just to MIP-start the model before returning it. The MIP-start is done by means of `MOI.set(model, VariablePrimalStart(), var_index, value)`."
		),
		Arg("allow-rotation", false,
			"The formulation will allow the pieces to rotate. Dummy rotated pieces are created and feeded to the discretization and model building. The only explicit difference in the formulation is that the demand is shared between the two rotations of a piece. Pieces that meet any of the following conditions do not have rotated counterparts: (i) the piece is a square; (ii) after rotation the piece does not fit into the 'original plate'/'large object'; (iii) the instance has another piece that is exactly the same as the rotated piece (other attributes like profit need to match too), in this case the demand of both (non-rotated) pieces is summed and shared."
		),
		Arg("mirror-plates", false,
			"Will error if the allow-rotation flag is not passed. If allow-rotation is enabled and mirror-plates is disabled, then the cut/plate enumeration is completely unaware of the fact (pieces are duplicated before enumeration and other changes are done at the formulation level). If both flags are enabled, then the enumeration is slightly changed: every plate can only have length smaller than or equal to the width (i.e., plates with length greater than width are rotated automatocally). This may cut by half the number of plates/constraints and, consequently, also reduce the number of cuts/variables (in other words, the whole size of the model). The Redundant-Cut improvement is disabled if this flag has any effect (because we lose the capability to distinguish between vertical and horizontal cuts)."
		),
		Arg("aggressive-hybridization", false,
			"Only can be used if hybridization-with-restricted is enabled. The default behavior of hybridization is to not increase the model size significantly (a small increase equivalent to the number of pieces may occur). To guarantee this, however, we need to abstain from hybridizations in a specific circumnstance. The circumnstance is: two or more pieces share the same length (or width), and the associated discretization point cannot be reached by a sum of smaller pieces lengths (or widths). Without increasing the number of variables (i.e., creating a different hybridized cut for each piece), we cannot hybridize for any of the pieces, because hybridizing for a single arbitrary piece would make the formulation lose the guarantee of optimality. Consequently, the default behaviour (i.e., without this flag) loses some opportunities for hybridization to keep optimality and to avoid considerable model size increase. This flag allows the formulation to seize these opportunities by creating multiple hybridized cuts, it also keeps the optimality but it can increase the model size considerably."
		)
	]
end

function Utilities.Args.throw_if_incompatible_options(::Val{:PPG2KP}, p_args)
	p_args["mirror-plates"] && !p_args["allow-rotation"] && throw(
		ArgumentError, "Flag mirror-plates cannot be used if allow-rotation" *
			" is disabled."
	)
	if p_args["aggressive-hybridization"] && !p_args["hybridize-with-restricted"]
		throw(
			ArgumentError, "Flag aggressive-hybridization cannot be used if" *
				"hybridize-with-restricted is disabled."
		)
	end
	alpha = p_args["pricing-alpha"]
	beta = p_args["pricing-beta"]
	alpha > 0.0 || alpha <= 1.0 || throw(ArgumentError(
		"Option pricing-alpha should follow: 0 < alpha <= 1.0, but is $(alpha)."
	))
	beta >= 0.0 || throw(ArgumentError(
		"Option pricing-beta must not be negative, but it is $(beta)."
	))
	beta <= 1.0 || @warn(
		"There is little sense in having a value above one in option" *
		" pricing-beta (it was $(beta)), are you sure you intended this?"
	)
	name = "pricing"
	Utilities.throw_if_unrecognized(
		name, p_args[name], ["expected", "furini", "becker", "none"]
	)
	name = "MIP-start"
	Utilities.throw_if_unrecognized(
		name, p_args[name], ["expected", "guaranteed", "none"]
	)
	building_time_limit = p_args["building-time-limit"]
	building_time_limit > 0.0 || throw(ArgumentError(
		"Option building-time-limit must be positive, but it is" *
		" $(building_time_limit)."
	))
	return
end

end # module

