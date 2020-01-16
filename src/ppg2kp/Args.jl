module Args

function get_arg_parse_settings()
	s = ArgParseSettings()
	ArgParse.add_arg_group(s, "PPG2KP-Specific Flags", "ppg2kp-specific-flags")
	@add_arg_table s begin
		"--ignore-d"
			help = "ignore the demand information during discretization, used to measure impact"
			nargs = 0
		"--use-c25"
			help = "add the tightening constraints 2.5"
			nargs = 0
		"--round2disc"
			help = "round the second child size to a discretized position"
			nargs = 0
		"--ignore-2th-dim"
			help = "ignore the dimension not being discretized during discretization, used to measure impact (does not affect flow)"
			nargs = 0
		"--final-pricing"
			help = "uses the best lb available (from --lower-bounds or --warm-start) to remove variables after solving the continuous relaxation"
			nargs = 0
		"--faithful2furini2016"
			help = "tries to be the most faithful possible to the description on the Furini2016 paper (more specifically the complete model, with the reductions, FOR NOW without a heuristic solution first and without pricing); the flags --no-cut-position, --no-redundant-cut and --no-furini-symmbreak disable parts of this reimplementation"
			nargs = 0
		"--no-redundant-cut"
			help = "disables Furini2016 Redundant-Cut reduction: a bunch of flags is used to check if some trim cut is necessary for optimality, or is dominated by other trim cuts"
			nargs = 0
		"--no-furini-symmbreak"
			help = "disables Furini2016 symmetry breaking: as trim cuts are needed in faithful2furini2016 mode, the symmetry-breaking is less restrictive than 'do not cut after midplate', it just removes each cut y > midplate in which plate_size - y is an existing cut; enabling this flag makes the code just use all discretized positions as cuts"
			nargs = 0
		"--no-cut-position"
			help = "disables Furini2016 Cut-Position reduction: if a trivial heuristic shows that no combination of six or more pieces fit a plate, then the plate may be cut with restricted cuts without loss to optimality"
			nargs = 0
		"--warm-start"
			help = "for now conflict with --faithful2furini2016) uses the heuristic described in Furini2016 to generate a initial primal feasible solution, warm-start the model with unrestricted cuts fixed to zero, and then unfix the unrestricted cuts to solve the complete model"
			nargs = 0
	end
	set_default_arg_group(s)
	s
end

end # module
