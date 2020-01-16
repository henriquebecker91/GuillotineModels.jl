module Args

using ArgParse

function parse_settings()
	s = ArgParseSettings()
	ArgParse.add_arg_group(s, "PPG2KP-Specific Flags", "ppg2kp-specific-flags")
	@add_arg_table s begin
		"--ignore-d"
			help = "Ignore the demand information during discretization, used to measure impact."
			nargs = 0
		"--use-c25"
			help = "Add the tightening constraints 2.5. Increase considerably the number of constraints hoping to improve the relaxation a little."
			nargs = 0
		"--round2disc"
			help = "Round the second child size to a discretized position."
			nargs = 0
		"--ignore-2th-dim"
			help = "Ignore the dimension not being discretized during discretization, used to measure impact."
			nargs = 0
		"--final-pricing"
			help = "Uses the best lb available (either from --lower-bounds or --warm-start) to remove variables after solving the continuous relaxation."
			nargs = 0
		"--faithful2furini2016"
			help = "Tries to be the most faithful possible to the description on the Furini2016 paper (more specifically the complete model, with the reductions, FOR NOW without a heuristic solution first and without pricing); the flags --no-cut-position, --no-redundant-cut and --no-furini-symmbreak disable parts of this reimplementation."
			nargs = 0
		"--no-redundant-cut"
			help = "Disables Furini2016 Redundant-Cut reduction: a bunch of flags is used to check if some trim cut is necessary for optimality, or is dominated by other trim cuts."
			nargs = 0
		"--no-furini-symmbreak"
			help = "Disables Furini2016 symmetry breaking: as trim cuts are needed in faithful2furini2016 mode, the symmetry-breaking is less restrictive than 'do not cut after midplate', it just removes each cut y > midplate in which plate_size - y is an existing cut; enabling this flag makes the code just use all discretized positions as cuts."
			nargs = 0
		"--no-cut-position"
			help = "Disables Furini2016 Cut-Position reduction: if a trivial heuristic shows that no combination of six or more pieces fit a plate, then the plate may be cut with restricted cuts without loss to optimality."
			nargs = 0
		"--warm-start"
			help = "(For now conflict with --faithful2furini2016) Uses the heuristic described in Furini2016 to generate a initial primal feasible solution, warm-start the model with unrestricted cuts fixed to zero, and then unfix the unrestricted cuts to solve the complete model."
			nargs = 0
	end
	set_default_arg_group(s)
	s
end

function check_flag_conflicts(p_args)
	@assert p_args["model"] == "PPG2KP"
	is_revised_furini = !p_args["faithful2furini2016"]
	p_args["final-pricing"] && !is_revised_furini && @error(
		"the final pricing technique is implemented just for Enhanced Furini" *
		" model as of now"
	)
	p_args["final-pricing"] && isempty(p_args["lower-bounds"]) &&
		!p_args["warm-start"] && @error(
		"the flag --final-pricing only makes sense if a lower bound is provided (either directly by --lower-bound or indirectly by --warm-start)"
	)
	p_args["final-pricing"] && p_args["relax2lp"] && @error(
		"the flags --final-pricing and --relax2lp should not be used together; it is not clear what they should do, and the best interpretation (solving the relaxed model and doing the final pricing, without solving the unrelaxed reduced model after) is not specially useful and need extra code to work that is not worth it"
	)
end
end # module
