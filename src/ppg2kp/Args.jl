module Args

using ArgParse
using ..Utilities.Args

function owned_args()
	return [
		Arg(
			"lower-bound", 0.0,
			"The use of this value is dependent on the other options selected, --final-pricing may use the values if --warm-start do not give better values."
		),
		Arg(
			"upper-bounds", Inf,
			"The use of these values is dependent on the other options selected, check code."
		),
		Arg(
			"round2disc", false,
			"Round the second child size to a discretized position."
		),
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
			"no-cut-position", false,
			"Disables Furini2016 Cut-Position reduction: if a trivial heuristic shows that no combination of six or more pieces fit a plate, then the plate may be cut with restricted cuts without loss to optimality."
		),
		Arg(
			"no-furini-symmbreak", false,
			"Disables Furini2016 symmetry breaking: as trim cuts are needed in faithful2furini2016 mode, the symmetry-breaking is less restrictive than 'do not cut after midplate', it just removes each cut y > midplate in which plate_size - y is an existing cut; enabling this flag makes the code just use all discretized positions as cuts."
		),
		Arg(
			"final-pricing", false,
			"Uses the best lb available (either from --lower-bounds or --warm-start) to remove variables after solving the continuous relaxation."
		),
		Arg(
			"warm-start", false,
			"(For now conflict with --faithful2furini2016) Uses the heuristic described in Furini2016 to generate a initial primal feasible solution, warm-start the model with unrestricted cuts fixed to zero, and then unfix the unrestricted cuts to solve the complete model."
		),
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
		)
	]
end

function accepted_args()
	return owned_args()
end

function parse_settings()
	s = ArgParseSettings()
	ArgParse.add_arg_group(s, "PPG2KP-Specific Flags", "ppg2kp-specific-flags")
	ArgParse.add_arg_table.(s, owned_args())
	#=
	@add_arg_table s begin
		"--lower-bounds"
			help = "Takes a single string, the string has N comma-separated-numbers where N is the number of instances passed. The use of these values is dependent on the other options selected, --final-pricing may use the values if --warm-start do not give better values."
			action = :store_arg
			arg_type = String
		"--upper-bounds"
			help = "Takes a single string, the string has N comma-separated-numbers where N is the number of instances passed, no spaces allowed. The use of these values is dependent on the other options selected, check code."
			action = :store_arg
			arg_type = String
		"--round2disc"
			help = "Round the second child size to a discretized position."
			action = :store_true
		"--faithful2furini2016"
			help = "Tries to be the most faithful possible to the description on the Furini2016 paper (more specifically the complete model, with the reductions, FOR NOW without a heuristic solution first and without pricing); the flags --no-cut-position, --no-redundant-cut and --no-furini-symmbreak disable parts of this reimplementation."
			action = :store_true
		"--no-redundant-cut"
			help = "Disables Furini2016 Redundant-Cut reduction: a bunch of flags is used to check if some trim cut is necessary for optimality, or is dominated by other trim cuts."
			action = :store_true
		"--no-cut-position"
			help = "Disables Furini2016 Cut-Position reduction: if a trivial heuristic shows that no combination of six or more pieces fit a plate, then the plate may be cut with restricted cuts without loss to optimality."
			action = :store_true
		"--no-furini-symmbreak"
			help = "Disables Furini2016 symmetry breaking: as trim cuts are needed in faithful2furini2016 mode, the symmetry-breaking is less restrictive than 'do not cut after midplate', it just removes each cut y > midplate in which plate_size - y is an existing cut; enabling this flag makes the code just use all discretized positions as cuts."
			action = :store_true
		"--final-pricing"
			help = "Uses the best lb available (either from --lower-bounds or --warm-start) to remove variables after solving the continuous relaxation."
			action = :store_true
		"--warm-start"
			help = "(For now conflict with --faithful2furini2016) Uses the heuristic described in Furini2016 to generate a initial primal feasible solution, warm-start the model with unrestricted cuts fixed to zero, and then unfix the unrestricted cuts to solve the complete model."
			action = :store_true
		"--ignore-d"
			help = "Ignore the demand information during discretization, used to measure impact."
			action = :store_true
		"--use-c25"
			help = "Add the tightening constraints 2.5. Increase considerably the number of constraints hoping to improve the relaxation a little."
			action = :store_true
		"--ignore-2th-dim"
			help = "Ignore the dimension not being discretized during discretization, used to measure impact."
			action = :store_true
	end
	=#
	set_default_arg_group(s)
	s
end

function check_flag_conflicts(p_args)
	is_revised_furini = !p_args["faithful2furini2016"]
	p_args["final-pricing"] && !is_revised_furini && @error(
		"the final pricing technique is implemented just for Enhanced Furini" *
		" model as of now"
	)
	p_args["final-pricing"] && isempty(p_args["lower-bounds"].vector) &&
		!p_args["warm-start"] && @error(
		"the flag --final-pricing only makes sense if a lower bound is provided (either directly by --lower-bound or indirectly by --warm-start)"
	)
end
end # module
