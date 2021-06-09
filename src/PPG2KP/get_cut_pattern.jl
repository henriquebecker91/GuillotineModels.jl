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
		if !iszero(fc)
			fc_idx = _consume_cut!(cuts_available, nz_cuts, fc)
			fc_idx !== nothing && push!(cut_idx_stack, fc_idx)
		end
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
# NOTE: only patterns and nzcs_vals are modified.
# Starting from the cuts closest to the piece extractions (i.e.,
# non-leaf nodes closest to a leaf node) start building the CutPattern
# (tree) bottom-up.
# The hybridizations are dealt with here.
function _bottom_up_tree_build!(
	patterns :: Dict{Int64, Vector{CutPattern{D, S}}},
	nz_cut_idx_stack,
	nz_cuts :: Vector{NTuple{3, P}},
	nz_cuts_ori :: BitArray{1},
	nz_cuts_pe :: Vector{D}, # The piece extractions from non-zero cuts
	nzcs_idxs :: Vector{P},
	nzcs_vals :: Vector{D},
	bmr :: ByproductPPG2KP{D, S, P},
	debug :: Bool = false
) :: Dict{Int64, Vector{CutPattern{D, S}}} where {D, S, P}
	for cut_idx in reverse(nz_cut_idx_stack)
		debug && @show cut_idx
		pp, fc, sc = nz_cuts[cut_idx]
		ce :: D = nz_cuts_pe[cut_idx]
		debug && @show pp, fc, sc, ce
		ppl, ppw = bmr.pli_lwb[pp][1], bmr.pli_lwb[pp][2]
		debug && @show ppl, ppw
		cut_ori = nz_cuts_ori[cut_idx]

		child_patts = Vector{CutPattern{D, S}}()
		# If `!iszero(ce)` then this cut is a double cut that extract piece
		# `ce` (only possible if flag `hybridize-with-restricted` is enabled).
		# However, if `isempty(nzcs_idx)` then no double cut ever
		# sold this piece type, the double cut was done just for trimming;
		# alternatively, if `iszero(nzcs_vals[only(nzcs_idx)])` then some
		# double cuts have sold pieces of type `ce` but these were already
		# processed, the remaining double cuts are to be considered as
		# trimming cuts. Finally, if the double cut really extracts a piece,
		# then we subtract it from `nzcs_vals`.
		#
		# We deal with double cuts in the following way: a dummy plate is
		# created to represent what would be the first child in a normal cut
		# but, in a double cut, it holds both the extracted piece (ce) and
		# the first child. If the cut is not a double cut, then we just add
		# the first child as the first child, instead of a dummy as the
		# first child.
		nzcs_idx = searchsorted(nzcs_idxs, ce)
		if isempty(nzcs_idx) || iszero(nzcs_vals[only(nzcs_idx)])
			ce = zero(D)
		else
			nzcs_vals[only(nzcs_idx)] -= one(D)
		end
		if !iszero(ce)
			# Get a safe size for the dummy: if there is not a first child
			# (it was killed by hubridization) then a safe size is just the
			# piece size; otherwise, we get the minimum necessary to pack
			# both the first child and the extracted piece.
			dummy_l, dummy_w = if iszero(fc)
				bmr.l[ce], bmr.w[ce]
			else
				if !cut_ori # if the cut inside the dummy is vertical
					max(bmr.l[ce], bmr.pli_lwb[fc][1]), bmr.w[ce] + bmr.pli_lwb[fc][2]
				else
					bmr.l[ce] + bmr.pli_lwb[fc][1], max(bmr.w[ce], bmr.pli_lwb[fc][2])
				end
			end
			dummy = CutPattern(
				# The cut orientation inside the dummy is the opposite as the
				# external cut (the first of the double cut): `!cut_ori` is the
				# orientation of the second (inner) cut (i.e., the cut in the dummy).
				dummy_l, dummy_w, !cut_ori,
				CutPattern{D, S}[CutPattern(bmr.l[ce], bmr.w[ce], ce)]
			)
			# Now just add the first child inside the dummy if it exists,
			# as we would do outside the dummy if there was no hybridization.
			if !iszero(fc) && haskey(patterns, fc) && !isempty(patterns[fc])
				push!(dummy.subpatterns, pop!(patterns[fc]))
				isempty(patterns[fc]) && delete!(patterns, fc)
			end
			# Finally, push the dummy as the first child of the cut/plate this
			# loop actually is about.
			push!(child_patts, dummy)
		else # If there is not double cut extraction just add the first child.
			# The fc can be zero yet because the current cut may be a double cut
			# that kills the first child and throws away the extracted piece.
			if !iszero(fc) && haskey(patterns, fc) && !isempty(patterns[fc])
				push!(child_patts, pop!(patterns[fc]))
				isempty(patterns[fc]) && delete!(patterns, fc)
			end
		end
		if !iszero(sc) && haskey(patterns, sc) && !isempty(patterns[sc])
			push!(child_patts, pop!(patterns[sc]))
			isempty(patterns[sc]) && delete!(patterns, sc)
		end
		!haskey(patterns, pp) && (patterns[pp] = Vector{CutPattern{D, S}}())
		push!(patterns[pp], CutPattern(
			ppl, ppw, cut_ori, child_patts
		))
	end

	patterns
end

# Internal function that does the heavy lifting with the data already
# extracted from the model.
function _get_cut_pattern(
	nzpe_idxs, nzpe_vals, nzcm_idxs, nzcm_vals, nzcs_idxs, nzcs_vals,
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
	# NOTE: if there are cut sells (from hybridize-with-restricted)
	# we cannot bail early, as the cuts themselves may extract/sell pieces.
	if isempty(nzpe_idxs) && isempty(nzcs_idxs)
		@assert isempty(nzpe_vals) && isempty(nzcs_vals)
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

	# Here we build the tree. All cuts are added, even if they do not
	# lead to extractions. If hybridize-with-restricted is enabled the
	# cuts themselves may extract pieces, so all of them will need to
	# be considered anyway. Before the call to _bottom_up_tree_build,
	# _nzcs_vals are all non-zero (as the 'nz' in the name says), after
	# the call the values must all be zero because the pieces sold from
	# double cuts extractions were all assigned to somewhere in the tree.
	# NOTE: sel_cuts is the subset of cuts given by nzcm_idxs, and the
	# cut_idx_stack was built using it so their indexes are not from the
	# universe set (as the values of nzcm_idxs) but from the subset
	# (the keys of sel_cuts and nzcm_idxs); sel_cuts_pe provides a way
	# to querying the piece extractions using these subset keys
	# (as bmr.cut_extraction understands thee universal keys).
	_nzcs_vals = copy(nzcs_vals)
	if isempty(nzcs_idxs)
		sel_cuts_pe = zeros(D, length(nzcm_idxs))
	else
		sel_cuts_pe = bmr.cut_extraction[nzcm_idxs]
	end
	@assert all(!iszero, _nzcs_vals)
	_bottom_up_tree_build!(
		patterns, cut_idx_stack, sel_cuts, ori_cuts, sel_cuts_pe,
		nzcs_idxs, _nzcs_vals, bmr, debug
	)
	@assert all(iszero, _nzcs_vals)

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
	rai2opi :: Vector{Union{D, Tuple{D, D}}},
	remaining_original_demand :: Vector{D}
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

		if iszero(remaining_original_demand[opi])
			@warn "The solution has more copies of the piece type number $opi" *
				" than exists demand for it. This is possible in G2OPP/G2CSP" *
				" because we treat demand as a lower bound for the solution to be" *
				" valid but we do not enforce equality."
		else
			remaining_original_demand[opi] -= 1
		end
		patterns[i] = CutPattern(e.length, e.width, opi)
	end

	return nothing
end

function _cp_rai2opi!(
	patterns :: Vector{CutPattern{D, S}},
	rad :: RotationAwareData{D, S}
) :: Nothing where {D, S}
	rai2opi = rad.sdi2opi[rad.rai2sdi]
	d = deepcopy(rad.original_demand)
	_inner_cp_rai2opi!(patterns, rai2opi, d)

	return nothing
end

# NOTE: exponential time on the number of children, but works over the
# non-simplified patterns that only have 2~3 children (2 if non-hybridized,
# possibly 3 if hybridized?).
# MAX and SUM are L and W if cuts_are_vertical and switched otherwise.
# max_s and sum_s are the length and width of all children plates/pieces
# if cuts_are_vertical and switched otherwise.
function _idxs_children_to_mirror_back_multiple(
	::Type{D}, MAX :: S, SUM :: S, max_s :: Vector{S}, sum_s :: Vector{S}
) :: Vector{D} where {D, S}
	@assert length(max_s) == length(sum_s)
	n = length(max_s)
	n > 8 && warn(
		"The code for unmirroring the solution is exponential on the number" *
		" immediate children of a pattern (i.e., the length of the subpatterns" *
		" field inside a CutPattern object). Found a subpattern with $n" *
		" children inside, this may slow down the computation."
	)
	for i in 1:(2^n - 1)
		turned = '1' .== collect(last(bitstring(i), n))

		new_max_s = ifelse.(turned, sum_s, max_s)
		new_sum_s = ifelse.(turned, max_s, sum_s)

		if maximum(new_max_s) <= MAX && sum(new_sum_s) <= SUM
			return (one(D):convert(D, n))[turned]
		end
	end

	@error "Tried to unmirror a subpattern but no combination of rotations" *
		" resulted in a valid solution. Orthogonal/max: $MAX and $max_s." *
		" Parallel/sum: $SUM and $sum_s."

	return D[]
end

function _idxs_children_to_mirror_back(
	children :: Vector{CutPattern{D, S}},
	L :: S, W :: S, cuts_are_vertical :: Bool
) :: Vector{D} where {D, S}
	isempty(children) && return D[]
	if isone(length(children))
		only_child = first(children)
		l, w = only_child.length, only_child.width
		if l > L || w > W # Some dimension does not fit, needs mirroring/rotation.
			if l > W && w > L # Impossible to fit the other way. Should not happen.
				@error "Could not un-mirror a single plate inside another." *
					"Plate $(L)x$(W) and inner plate/piece $(l)x$(w)."
				return D[]
			else
				return D[one(D)]
			end
		else # The only child fits, no need to rotate.
			return D[]
		end
	else # length(patterns) > 1
		cw = getfield.(children, :width)
		cl = getfield.(children, :length)
		if cuts_are_vertical
			sum_width = sum(cw)
			max_length = maximum(cl)
			if sum_width > W || max_length > L
				return _idxs_children_to_mirror_back_multiple(D, L, W, cl, cw)
			else # everything fits
				return D[]
			end
		else # !cuts_are_vertical
			sum_length = sum(cl)
			max_width = maximum(cw)
			if sum_length > L || max_width > W
				return _idxs_children_to_mirror_back_multiple(D, W, L, cw, cl)
			else # everything fits
				return D[]
			end
		end
	end
end

function _rotate_pattern(
	p :: CutPattern{D, S}
) :: CutPattern{D, S} where {D, S}
	return CutPattern(
		p.width, p.length, p.piece_idx,
		isempty(p.subpatterns) ? p.cuts_are_vertical : !p.cuts_are_vertical,
		p.subpatterns # we do only rotate the subpatterns if necessary, so not here
	)
end

function _rec_try_mirror_back!(
	p :: CutPattern{D, S}
) :: Nothing where {D, S}
	# Get the index of the subpatterns that will be rotated.
	c_idxs = _idxs_children_to_mirror_back(
		p.subpatterns, p.length, p.width, p.cuts_are_vertical
	)
	# Rotate the subpatterns (but not their subpatterns).
	p.subpatterns[c_idxs] .= _rotate_pattern.(p.subpatterns[c_idxs])
	# Now recurse for each subpattern, rotated or not.
	foreach(_rec_try_mirror_back!, p.subpatterns)

	return
end

# A new CutPattern is returned and the old one may be invalid and share
# memory with the new one (so nothing should use it anymore).
function _try_mirror_back!(
	pattern :: CutPattern{D, S}, L :: S, W :: S
) :: CutPattern{D, S} where {D, S}
	L_, W_ = pattern.length, pattern.width
	if L_ > L || W_ > W
		if pattern.length > W || pattern.width > L
			@warn "While unmirroring the solution failed to fit the root pattern." *
				"Original dimensions? $L and $W. Root pattern? $L_ and $W_."
		end
		pattern = _rotate_pattern(pattern)
	end

	_rec_try_mirror_back!(pattern)

	return pattern
end

function _single_pattern(
	patterns :: AbstractVector{CutPattern{D, S}}
) :: Union{CutPattern{D, S}, Nothing} where {D, S}
	qt_patterns = length(patterns)
	if qt_patterns > 1
		error(
			"The problem only admits a single pattern but multiple were found." *
			" The patterns follow:\n" * join(to_pretty_str.(patterns), "\n") * "\n"
		)
	end
	return iszero(qt_patterns) ? nothing : only(patterns)
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

	# TODO: the fact dc_sells only exist in hybridized models should probably
	# be regarded as a implementation detail, and we should use a boolean
	# in the bmr.preprocess_byproduct instead. But, for now, this will
	# suffice.
	hybridize_with_restricted = haskey(JuMP.object_dictionary(model), :dc_sells)
	nzcs_idxs, nzcs_vals = if hybridize_with_restricted
		cs = model[:dc_sells] # Cut Sells
		gather_nonzero(cs, D)
	else
		D[], D[]
	end

	if debug
		println("Start of: formulation info before _get_cut_pattern.")
		@show nzpe_idxs
		@show nzpe_vals
		@show nzcm_idxs
		@show nzcm_vals
		if hybridize_with_restricted
			@show bmr.preprocess_byproduct.cut_extraction[nzcm_idxs]
			@show value.(cs[nzcs_idxs])
		end
		@show nzcs_idxs
		@show nzcs_vals
		@show value.(pe[nzpe_idxs])
		@show value.(cm[nzcm_idxs])
		println("End of: formulation info before _get_cut_pattern.")
	end

	# Call the method that deals only with the data, and not with the JuMP.Model.
	# NOTE: this can alter the arrays above.
	patterns = _get_cut_pattern(
		nzpe_idxs, nzpe_vals, nzcm_idxs, nzcm_vals, nzcs_idxs, nzcs_vals,
		bmr.preprocess_byproduct, debug
	) :: Vector{CutPattern{D, S}}

	if debug
		println("Start of: formulation info after _get_cut_pattern.")
		@show nzpe_idxs
		@show nzpe_vals
		@show nzcm_idxs
		@show nzcm_vals
		if hybridize_with_restricted
			@show bmr.preprocess_byproduct.cut_extraction[nzcm_idxs]
			@show value.(cs[nzcs_idxs])
		end
		@show nzcs_idxs
		@show nzcs_vals
		@show value.(pe[nzpe_idxs])
		@show value.(cm[nzcm_idxs])
		println("End of: formulation info after _get_cut_pattern.")
	end

	# All the CutPattern extraction procedure is rotation-unaware. If rotation
	# is allowed, the piece indexes in the leaf nodes of CutPattern need to be
	# changed to the original piece indexes.
	if bmr.rad !== nothing
		_cp_rai2opi!(patterns, bmr.rad :: RotationAwareData{D, S})
		if bmr.preprocess_byproduct.mirror_plates
			L, W = bmr.preprocess_byproduct.L, bmr.preprocess_byproduct.W
			patterns = _try_mirror_back!.(patterns, L, W)
		end
	end

	is_single_pattern = isa(problem, Val{:G2KP}) || isa(problem, Val{:G2OPP})
	if is_single_pattern
		return _single_pattern(patterns) :: Union{CutPattern{D, S}, Nothing}
	end
	return patterns :: Vector{CutPattern{D, S}}
end
