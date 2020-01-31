module KnapsackPlates

module Args
	import ...Utilities
	import ...Utilities.Args: Arg
	#export accepted_arg_list, throw_if_incompatible_options
	# No extra flags for this model yet.
	Utilities.Args.accepted_arg_list(::Val{:KnapsackPlates}) = Vector{Arg}()
	function Utilities.Args.throw_if_incompatible_options(
		::Val{:KnapsackPlates}, p_args
	)
		nothing
	end
end

using JuMP

import ..build_model
"""
    build_model(::Val{:KnapsackPlates}, model, d, p, l, w, L, W, [options])

The KnapsackPlates models each plate as a 1D knapsack of pieces and smaller
subplates. Basically, a plate is allowed to be cut in just one orientation
(except the original plate, that has a hack to allow the first stage to
be any of the two orientations) and, this way, the plate is treated as a
1D knapsack problem using the size of pieces and plates in the considered
dimension.

NOTE: THIS MODEL IS BROKEN. The source has the skeleton of the implementation
but it fails compilation as it is.

```
TODO: add high-level description of the mathematical model here.
```

"""
function build_model(
	::Val{:KnapsackPlates}, model, d :: Vector{D}, p :: Vector{P},
	l :: Vector{S}, w :: Vector{S}, L :: S, W :: S,
	options :: Dict{String, Any} = Dict{String, Any}()
) where {D, S, P}
	@assert length(d) == length(l) && length(l) == length(w)
	num_fplate_types = convert(D, length(d))

	@variable(model, 0 <= fplate[i = 1:num_fplate_types] <= d[i], Int)
	@variable(model, 0 <= iplate[i = 1:num_iplate_types] <= qt_fit_op[i], Int)

	@objective(model, Max,
		sum(p[fpi] * fplate[fpi] for fpi = 1:num_fplate_types)
	)

	# The original plate issue (no way to cut allow cutting both horizonatlly and
	# vertically) may be solved by having a default cut orientation for the
	# original plate, and allowing to cut a same-size-but-opposite-orientation
	# plate. Such plate has length L and is contained in the list of plates to
	# be cut vertically.
	@constraint(model,
		sum(lfp .* fp) + sum(lip2cv .* ip2cv) <= L
	)
	# iip2ch: Index of Intermediary Plates To Cut Horizontally
	for iip2ch = 1:num_ip2ch_types
		fitting_fp = fp_fitting_iipch[iip2ch]
		fitting_ip2cv = iip2cv_fitting_iip2ch[iip2ch]
		isempty(fitting_fp) || isempty(fitting_ip2cv) && continue
		@constraint(model,
			sum(lfp[fitting_fp] .* fp[fitting_fp]) +
			sum(lip2cv[fitting_ip2cv] .* ip2cv[fitting_ip2cv]) <=
			lip2ch[iip2ch] * ip2ch[iip2ch]
		)
	end

	for iip2cv = 1:num_ip2cv_types
		fitting_fp = fp_fitting_iipcv[iip2cv]
		fitting_ip2ch = iip2ch_fitting_iip2cv[iip2cv]
		isempty(fitting_fp) || isempty(fitting_ip2ch) && continue
		@constraint(model,
			sum(wfp[fitting_fp] .* fp[fitting_fp]) +
			sum(wip2ch[fitting_ip2ch] .* ip2ch[fitting_ip2ch]) <=
			wip2cv[iip2cv] * ip2cv[iip2cv]
		)
	end

	model
end

end # module
