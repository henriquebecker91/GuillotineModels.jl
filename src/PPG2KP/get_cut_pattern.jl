# INTERNAL METHOD USED ONLY IN get_cut_pattern
# If the "pattern" is the extraction of a single piece return either:
# (1) a CutPattern representing the piece, if it is the same size as the
# original plate; otherwise (2) a CutPattern of the original plate containing
# a CutPattern representing the piece.
function _extraction_pattern(
	bmr :: ByproductPPG2KP{D, S, P}, e_idx
) :: CutPattern{D, S} where {D, S, P}
	pli, pii = bmr.np[e_idx] # the plate index and the piece index
	L, W = bmr.pli_lwb[pli] # the plate dimensions
	piece = CutPattern(bmr.l[pii], bmr.w[pii], pii)
	L == bmr.l[pii] && W == bmr.w[pii] && return piece
	return CutPattern(L, W, false, CutPattern{D, S}[piece])
end

# INTERNAL USE.
# For each large object that has a single piece extracted directly from it,
# i.e., no intermediary plates or cuts, return the indexes in `np` of such
# extractions; also, mutate nzpe_idxs and nzcm_vals to remove the positions
# refering to extractions directly from the large object.
function _single_so_from_lo_extractions!(
	np :: Vector{Tuple{P, D}}, nzpe_idxs :: Vector{P}, nzpe_vals :: Vector{D}
) :: Vector{P} where {D, P}
	to_delete = D[]
	extractions = P[]
	for (nzpe_idx, np_idx) in pairs(nzpe_idxs)
		pli = np[np_idx][1]
		!isone(pli) && continue
		for _ = 1:nzpe_vals[nzpe_idx]
			push!(extractions, np_idx)
		end
		push!(to_delete, nzpe_idx)
	end
	deleteat!(nzpe_idxs, to_delete)
	deleteat!(nzpe_vals, to_delete)

	return extractions
end

# INTERNAL USE.
# Find the index in `cuts` of the first cut that has `pp` as its parent plate
# and has a positive `num_uses_in_sol` too; return the index after
# decrementing the `num_uses_in_sol[index]` by one. Returns `nothing` if
# there is no plate that satisfy such conditions.
# Note: `cuts` is not changed, only `num_uses_in_sol` is changed
function _consume_cut!(
   num_uses_in_sol :: Vector{D}, cuts :: Vector{NTuple{3, P}}, pp :: P
) :: Union{Nothing, P} where {D, P}
	for i in keys(cuts)
		if num_uses_in_sol[i] > zero(D) && first(cuts[i]) == pp
			num_uses_in_sol[i] -= one(D)
			return convert(P, i)
		end
	end
	return nothing
end

# INTERNAL USE.
# Given a list of the cuts used in a solution, the number of times their appear
# in the solution, and the index of the root cut in such list, return a list of
# cut indexes in topological ordering (any non-root cut over some plate only
# appears if a previous cut has made a copy of that plate type available).
# SEE: https://discourse.julialang.org/t/unreachable-reached-at-0x7fa478093547-in-julia-1-0-5/36404
function _build_cut_idx_stack(
	nz_cuts :: Vector{NTuple{3, P}},
	qt_cuts :: Vector{D},
	root_cut_idxs :: Vector{P},
	debug :: Bool = false
) :: Vector{P} where {D, P}
	cut_idx_stack = Vector{P}()
	append!(cut_idx_stack, expand(qt_cuts[root_cut_idxs], root_cut_idxs))
	cuts_available = deepcopy(qt_cuts)
	cuts_available[root_cut_idxs] .= zero(D)

	next_cut = one(P)
	while next_cut <= length(cut_idx_stack)
		_, fc, sc = nz_cuts[cut_idx_stack[next_cut]]
		@assert !iszero(fc)
		fc_idx = _consume_cut!(cuts_available, nz_cuts, fc)
		fc_idx !== nothing && push!(cut_idx_stack, fc_idx)
		if !iszero(sc)
			sc_idx = _consume_cut!(cuts_available, nz_cuts, sc)
			sc_idx !== nothing && push!(cut_idx_stack, sc_idx)
		end
		next_cut += one(P)
	end
	# If cuts_available is different from a vector of zeros then some cuts
	# were not consumed (I am not sure if this should be possible).
	debug && @show cuts_available

	return cut_idx_stack
end

# INTERNAL USE.
# NOTE: only patterns is modified.
# Given the non-zero (i.e., used) piece extractions (and the number of
# times they are used), add them to the patterns dictionary (that
# translates a plate id to a list of patterns for which the root node
# is a plate of that type, in this case, we are adding the patterns
# that are just piece extractions).
function _add_used_extractions!(
	patterns :: Dict{Int64, Vector{CutPattern{D, S}}},
	nzpe_idxs, nzpe_vals, bmr :: ByproductPPG2KP{D, S, P},
	debug :: Bool = false
) :: Dict{Int64, Vector{CutPattern{D, S}}} where {D, S, P}
	for (i, np_idx) in pairs(nzpe_idxs)
		pli, pii = bmr.np[np_idx]
		if debug
			@show np_idx
			@show pli, pii
			@show nzpe_vals[i]
			@show bmr.l[pii], bmr.w[pii]
		end
		extractions = _extraction_pattern.(bmr, repeat([np_idx], nzpe_vals[i]))
		if haskey(patterns, pli)
			append!(patterns[pli], extractions)
		else
			patterns[pli] = extractions
		end
	end
	patterns
end

# INTERNAL USE.
# NOTE: only patterns is modified.
# Starting from the cuts closest to the piece extractions (i.e.,
# non-leaf nodes closest to a leaf node) start building the CutPattern
# (tree) bottom-up.
function _bottom_up_tree_build!(
	patterns :: Dict{Int64, Vector{CutPattern{D, S}}},
	nz_cut_idx_stack,
	nz_cuts :: Vector{NTuple{3, P}},
	nz_cuts_ori :: BitArray{1},
	bmr :: ByproductPPG2KP{D, S, P},
	debug :: Bool = false
) :: Dict{Int64, Vector{CutPattern{D, S}}} where {D, S, P}
	for cut_idx in reverse(nz_cut_idx_stack)
		debug && @show cut_idx
		pp, fc, sc = nz_cuts[cut_idx]
		debug && @show pp, fc, sc
		child_patts = Vector{CutPattern{D, S}}()
		if !iszero(fc) && haskey(patterns, fc) && !isempty(patterns[fc])
			push!(child_patts, pop!(patterns[fc]))
			isempty(patterns[fc]) && delete!(patterns, fc)
		end
		if !iszero(sc) && haskey(patterns, sc) && !isempty(patterns[sc])
			push!(child_patts, pop!(patterns[sc]))
			isempty(patterns[sc]) && delete!(patterns, sc)
		end
		ppl, ppw = bmr.pli_lwb[pp][1], bmr.pli_lwb[pp][2]
		debug && @show ppl, ppw
		!haskey(patterns, pp) && (patterns[pp] = Vector{CutPattern{D, S}}())
		push!(patterns[pp], CutPattern(
			ppl, ppw, nz_cuts_ori[cut_idx], child_patts
		))
	end

	patterns
end

# Internal function that does the heavy lifting with the data already
# extracted from the model.
function _get_cut_pattern(
	nzpe_idxs, nzpe_vals, nzcm_idxs, nzcm_vals,
	bmr :: ByproductPPG2KP{D, S}, debug :: Bool
) :: Vector{CutPattern{D, S}} where {D, S}
	# 1. If there are extractions of small objects from large objects,
	#    we keep them in a separate array and remove these extractions from
	#    the data (the rest of the code then can safely assume the large
	#    objects are either cut or unused, not directly extracted).
	# 2. If not all extractions were made directly from large objects, then the
	#    cuts made over large objects are in `cuts_made`, find them (i.e.,
	#    traverse non-zero cuts_made variables checking for a hcut or vcut that
	#    has the original plate as the parent plate).
	# 3. Create a topological ordering of the cuts starting from the root cuts.
	# 4. Initialize a pool of plate type "uses"/patterns with the piece
	#    extractions in the solution (i.e., leaf nodes).
	# 5. Traverse the reverse of the topological ordering and build the
	#    'pattern tree' in a bottom-up fashion. The children of each pattern
	#    are queryed from the pool of plate type "uses", and removed from it
	#    to be only present inside their parent pattern. At the end of the
	#    process there should remain only patterns starting from an original
	#    plate (i.e., large object) in the pool.
	lo2so_idxs = _single_so_from_lo_extractions!(bmr.np, nzpe_idxs, nzpe_vals)
	lo2so_patterns = _extraction_pattern.(
		(bmr,), lo2so_idxs
	) :: Vector{CutPattern{D, S}}
	@assert all(!isequal(zero(eltype(nzpe_idxs))), nzpe_idxs)
	@assert all(!isequal(zero(D)), nzpe_vals)
	# If the call to _single_so_from_lo_extractions has exhausted the
	# list of extractions, then there is nothing else to do (i.e.,
	# lo2so_patterns has the solution already).
	if isempty(nzpe_idxs)
		@assert isempty(nzpe_vals)
		return lo2so_patterns
	end

	# The cuts actually used in the solution.
	sel_cuts = bmr.cuts[nzcm_idxs]
	debug && @show sel_cuts
	# If the cut in `sel_cuts` is vertical or not.
	ori_cuts = nzcm_idxs .>= bmr.first_vertical_cut_idx
	# The index of the root cuts (cut over the original plates) in `sel_cuts`.
	root_idxs = findall(cut -> isone(cut[1]), sel_cuts)
	isempty(root_idxs) && return CutPattern{D, S}[]

	cut_idx_stack = _build_cut_idx_stack(sel_cuts, nzcm_vals, root_idxs, debug)

	# `patterns` translates a plate index (pli) into a list of all "uses" of that
	# plate type. Such "uses" are CutPattern objects, either pieces or more
	# complex patterns. When we discover a plate is used in the solution such
	# "use" is added to the known uses in `patterns`, and we arbitrarily
	# associate it to some "uses" of its children (which are then are removed
	# from `patterns`). Consequently, at the end, the `patterns` should have only
	# the key of the original plate, and its associated value should be a Vector
	# of a single CutPattern.
	patterns = Dict{Int64, Vector{CutPattern{D, S}}}()

	# Insert all piece extractions (i.e., plate to piece patterns) into
	# `patterns`.
	_add_used_extractions!(patterns, nzpe_idxs, nzpe_vals, bmr, debug)

	_bottom_up_tree_build!(patterns, cut_idx_stack, sel_cuts, ori_cuts, bmr, debug)

	if !isone(length(patterns)) || !haskey(patterns, 1)# || !isone(length(patterns[1]))
		println("BUG AT GET_CUT_PATTERN")
		for (key, subpatts) in patterns
			@show key
			for subpatt in subpatts
				println(to_pretty_str(subpatt))
			end
		end
	end

	@assert isone(length(patterns))
	@assert haskey(patterns, 1)
	#@assert isone(length(patterns[1]))

	return vcat(lo2so_patterns, patterns[1])
end

function _inner_cp_rai2opi!(
	patterns :: AbstractVector{CutPattern{D, S}},
	rai2opi :: AbstractVector{Union{D, Tuple{D, D}}},
	remaining_original_demand :: AbstractVector{D}
) :: Nothing where {D, S}
	for (i, e) in pairs(patterns)
		if iszero(e.piece_idx) # not a piece extraction
			_inner_cp_rai2opi!(e.subpatterns, rai2opi, remaining_original_demand)
			continue
		end

		opis = rai2opi[e.piece_idx] # original piece index(es)
		opi :: D = if isa(opis, D) # The dummy maps to a single original piece.
			opis
		else # The dummy maps to two original pieces (need to manage demand).
			iszero(remaining_original_demand[opis[1]]) ? opis[2] : opis[1]
		end

		@assert !iszero(remaining_original_demand[opi])
		remaining_original_demand[opi] -= 1
		patterns[i] = CutPattern(e.length, e.width, opi)
	end

	return nothing
end

function _cp_rai2opi!(
	patterns :: AbstractVector{CutPattern{D, S}},
	rad :: RotationAwareData{D, S}
) :: Nothing where {D, S}
	rai2opi = rad.sdi2opi[rad.rai2sdi]
	d = deepcopy(rad.original_demand)
	_inner_cp_rai2opi!(patterns, rai2opi, d)

	return nothing
end

import ..get_cut_pattern
#=
@timeit TIMER function get_cut_pattern(
	model_type :: Val{:PPG2KP}, model :: JuMP.Model, ::Type{D}, ::Type{S},
	build_model_return :: ModelByproduct{D, S, P}
) :: CutPattern{D, S} where {D, S, P}
=#
@timeit TIMER function get_cut_pattern(
	problem :: SIMILAR_4, formulation :: Val{:PPG2KP}, model :: JuMP.Model,
	build_model_return :: ModelByproduct{D, S, P}
) where {D, S, P}
	# local constant to alternate debug mode (method will not take a debug flag)
	debug = false

	bmr = build_model_return

	bmr.found_optimum && return bmr.optimum_if_found

	pe = model[:picuts] # Piece Extractions
	cm = model[:cuts_made] # Cuts Made

	# non-zero {piece extractions, cuts made} {indexes,values}
	nzpe_idxs, nzpe_vals = gather_nonzero(pe, D)
	nzcm_idxs, nzcm_vals = gather_nonzero(cm, D)
	if debug
		println("Start of: formulation info before _get_cut_pattern.")
		@show nzpe_idxs
		@show nzpe_vals
		@show nzcm_idxs
		@show nzcm_vals
		@show value.(pe[nzpe_idxs])
		@show value.(cm[nzcm_idxs])
		println("End of: formulation info before _get_cut_pattern.")
	end

	# NOTE: this can alter the arrays above.
	patterns = _get_cut_pattern(
		nzpe_idxs, nzpe_vals, nzcm_idxs, nzcm_vals, bmr.preprocess_byproduct, debug
	)

	if debug
		println("Start of: formulation info after _get_cut_pattern.")
		@show nzpe_idxs
		@show nzpe_vals
		@show nzcm_idxs
		@show nzcm_vals
		@show value.(pe[nzpe_idxs])
		@show value.(cm[nzcm_idxs])
		println("End of: formulation info after _get_cut_pattern.")
	end

	# All the CutPattern extraction procedure is rotation-unaware. If rotation
	# is allowed, the piece indexes in the leaf nodes of CutPattern need to be
	# changed to the original piece indexes.
	bmr.rad !== nothing && _cp_rai2opi!(patterns, bmr.rad)

	# Call the method that deals only with the data, and not with the JuMP.Model.
	is_single_pattern = isa(problem, Val{:G2KP}) || isa(problem, Val{:G2OPP})
	is_single_pattern && return only(patterns)
	return patterns
end
