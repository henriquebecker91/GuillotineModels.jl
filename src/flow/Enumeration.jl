module Enumeration

export Node, Edge, gen_all_edges

"A node of the GuillotineModels.Flow model."
struct Node{S, N}
	"Unique identifier of the node in the model."
	idx :: N
	"The size of the plate dimension parallel to the cuts."
	par :: S
	"The size of the plate dimension perpendicular to the cuts."
	per :: S
	"The cut orientation: one means vertical cuts; two means horizontal cuts."
	ori :: UInt8
end

"An edge of the GuillotineModels.Flow Model."
struct Edge{N, E}
  "Unique identifier of the edge in the model."
	indx :: E
	"The node from where the flow comes. Source."
	head :: N
	"The node from where the flow ends. Sink."
	tail :: N
	"The backward edge global index, if it exists, otherwise zero."
	back :: E
end

"""
    gen_rr_fow_edges!(::Type{N}, d :: [D], l :: [S], L :: S) :: ([N], [(S, S)])

# Arguments

An integer type capable of representing the nodes `N`. The pieces demand `d`,
their length (or width) `l`, and a plate length (or width) `L`. This method
assumes the selected pieces fit the selected plate (on both dimensions, this
is, in the one that was given and the one that wasn't).

# Return

A tuple of two vectors:
1) a vector `v` of length `L`; `v[y] == 1` iff there is a linear combination
  of the pieces (respecting demand) that gives `y`; otherwise `v[y] == 0`.
2) a vector containing all ordered pairs of linear combinations in which the
  difference between the two is exactly the size of a piece.
"""
function gen_rr_fow_edges!(
	::Type{N}, # eltype of the marked vector
	d :: Vector{D}, l :: Vector{S}, L :: S
) :: Tuple{Vector{N}, Vector{Tuple{S, S}}} where {D, S, N}
	@assert length(d) == length(l)
	n = length(d)
	marked :: Vector{N} = fill(zero(N), L)
	edges = Vector{Tuple{S, S}}() # start point and end point
	iszero(n) && return marked, edges
	# Mark the cuts of the first piece.
	marked[l[1]] = one(N)
	push!(edges, (zero(S), l[1]))
	y = l[1] # y: used to iterate capacity, inherited from knapsack papers
	for _ = 2:min(d[1], L ÷ l[1])
		y_ = y
		push!(edges, (y_, y += l[1]))
		marked[y] = one(N)
	end
	# Mark the cuts of all other pieces.
	for pii = 2:n # PIece Index
		li = l[pii] # length of i
		di = d[pii] # demand of i
		for y = (L - li):-1:1
			if !iszero(marked[y])
				yrli = y + li # yrli: y + repeated lengths of i
				push!(edges, (y, yrli))
				marked[yrli] = one(N)
				for _ = 2:di
					yrli_ = yrli
					yrli += li
					yrli > L && break
					push!(edges, (yrli_, yrli))
					marked[yrli] = one(N)
				end
			end
		end
		# Mark cuts using just multiple copies of the same piece type (could be
		# done using a dummy cut at index zero with value true, but the array has
		# no index zero).
		y = li
		push!(edges, (zero(S), y))
		marked[y] = one(N)
		for _ = 2:min(di, L ÷ li)
			y_ = y
			push!(edges, (y_, y += li))
			marked[y] = one(N)
		end
	end

	marked, edges
end

"""
    reduce_dlw(d::[D], l::[S], w::[S], L::S, W::S) :: ([D], [S], [S], D)

Create a copy of `d`, `l`, and `w` with only the pieces that fit the plate
`L` x `W` (i.e., pieces with both dimensions smaller than the corresponding
plate dimension). The last value returned is the number of pieces that pass
this criteria (and, consequently, the common size of the returned vectors).
"""
function reduce_dlw(
	d :: Vector{D}, l :: Vector{S}, w :: Vector{S}, L :: S, W :: S
) :: Tuple{Vector{D}, Vector{S}, Vector{S}, D} where {D, S}
	@assert length(d) == length(l)
	@assert length(d) == length(w)
	N = convert(D, length(d))
	# Remove items that cannot fit into L anyway.
	d_, l_, w_ = similar.((d, l, w))
	j = zero(D)
	for i in one(D):N
		if l[i] <= L && w[i] <= W
			j += one(D)
			d_[j] = d[i]
			l_[j] = l[i]
			w_[j] = w[i]
		end
	end
	resize!(d_, j)
	resize!(l_, j)
	resize!(w_, j)
	d_, l_, w_, j
end

"""
    merge_duplicates(d :: [D], per :: [S]) :: ([D], [S])

If `per` has no repeated values just return the parameters. Otherwise,
for each value in `per` with duplicates, keep just one copy of such value
inside `per` and delete the other copies, and the same positions in `d`;
the position in `d` corresponding to the kept copy will have the sum of
the `d` values for all copies deleted and the copy kept.
Sort `per` by increasing value and swap positions in `d` to keep the
correspondence.
"""
function merge_duplicates(
	d :: Vector{D}, per :: Vector{S}
) :: Tuple{Vector{D}, Vector{S}} where {D, S}
	@assert length(d) == length(per)
	isempty(d) && return d, per

	y2d = fill(zero(D), maximum(per))
	has_duplicates = false
	for i in eachindex(d)
		!iszero(y2d[per[i]]) && (has_duplicates = true)
		y2d[per[i]] += d[i]
	end

	!has_duplicates && return d, per

	per_ = Vector{S}()
	d_ = Vector{D}()
	for (y, dy) in enumerate(y2d)
		if !iszero(dy)
			push!(d_, dy)
			push!(per_, y)
		end
	end

	d_, per_
end

"""
    globalize!(...)

The model may be seen as constituted of many independent graphs, but as there
is interaction between them (i.e., subplate edges increase the flow of a
backward edge in other graph) their nodes and edges must follow a global
numbering and the backward edges must be known by all of graphs (as subplate
edges need to refer to them as child plates).

This method takes the "local/raw" lists of node (`y2node_idx`) and edge
(`raw_fow_edges`) numbering (that is just the normal cut positions); the index
of the last globalized node (`last_gnode_idx`) and edge (`last_gedge_idx`); a
matrix that gives the piece index given the piece dimensions (`lw2pii`); a
tridimensional array that gives the backward edge index given its plate
dimensions and allowed cut orientation (`ppo2gbedge_idx`); the allowed cut
`ori`entation of the given graph; and the size of the plate dimension parallel
to the allowed cuts (`PAR`).

The method return the "globalized" list of nodes and edges, as the index of
the last globalized node and edge.

The method changes `ppo2gbedge_idx`: if a subplate (forward) edge is created
and the corresponding backward edge was not yet numbered, it is given a
"globalized" index (all backward edges are created at another step at the end).
"""
function globalize!(
	y2node_idx     :: Vector{N},
	raw_fow_edges  :: Vector{Tuple{S, S}}, # raw forward edges
	last_gnode_idx :: N,
	last_gedge_idx :: E,
	lw2pii         :: AbstractArray{D, 2},
	ppo2gbedge_idx :: Array{E, 3},
	ori            :: UInt8,
	PAR            :: S
) :: Tuple{Vector{Node{S, N}}, Vector{Edge{N, E}}, N, E} where {D, S, N, E}
	# The first node of a flow is the point zero, that does not appear in a
	# discretization because the array has no position zero.
	lgni = origin = last_gnode_idx + one(N)
	lgei = last_gedge_idx

	glo_nodes = Vector{Node{S, N}}()
	push!(glo_nodes, Node{S, N}(origin, PAR, zero(S), ori))

	for y in one(S):length(y2node_idx)
		if !iszero(y2node_idx[y])
			y2node_idx[y] = (lgni += one(N))
			push!(glo_nodes, Node{S, N}(lgni, PAR, y, ori))
		end
	end

	# glo_fow_edges: vector in which the values are the globalized forward edges
	# (we call it globalized because indx is the global name of that arc, and
	# both head and tail are global names of nodes).
	glo_fow_edges = Vector{Edge{N, E}}()

	for (y1, y2) in raw_fow_edges
		@assert y1 >= zero(y1) # there are no negative positions
		@assert y2 > y1 # i.e., all edges are forward
		indx = (lgei += one(E))
		if iszero(y1)
			head = origin
		else
			head = y2node_idx[y1]
		end
		tail = y2node_idx[y2]
		# NOTE: the backward edges are not created here. The backward edges are
		# created after all forward edges in the whole model are created.
		# NOTE: the ppo2gbedge_idx is indexed by parallel-perpendicular-orientation
		# but as we check a child plate where the orientation is inversed, then
		# the ori value is inverted, and parallel and perpendicular for the current
		# flow/plate are reversed for the child plate.
		bedge_ppo = (y2 - y1, PAR, 0x03 - ori)
		if iszero(ppo2gbedge_idx[bedge_ppo...])
			ppo2gbedge_idx[bedge_ppo...] = (lgei += one(E))
		end
		back = ppo2gbedge_idx[bedge_ppo...]
		edge = Edge{N, E}(indx, head, tail, back)
		#(iszero(head) || iszero(tail)) && @show edge, @__LINE__
		push!(glo_fow_edges, edge)
		# If this foward edge cuts a plate with the exact same size as some piece
		# then we add the "final plate"/piece edge too.
		yd = y2 - y1
		max_par = size(lw2pii, 1)
		max_per = size(lw2pii, 2)
		#@show PAR, yd, max_par, max_per, ori
		if PAR <= max_par && yd <= max_per && !iszero(lw2pii[PAR, yd])
			indx = (lgei += one(E))
			edge = Edge{N, E}(indx, head, tail, lw2pii[PAR, yd])
			push!(glo_fow_edges, edge)
		end
	end

	glo_nodes, glo_fow_edges, lgni, lgei
end

"""
    gen_u_fow_edges(glo_nodes, y2node_idx, last_gedge_idx, ppo2gbedge_idx)

This method generate the unrestricted subplate edges, i.e., the subplate edges
representing a cut that creates a subplate in which the size of the dimension
that is perpendicular to the cut is not the size of a single piece, but of a
linear combination of them. Such edges are necessary to solve the unrestricted
case but not to solve restricted case.

# Arguments

* `glo_nodes::Vector{Node{S, N}}`: A list of globalized nodes from a
  single graph (asserts will fail if nodes are from different graphs).
* `y2node_idx::Vector{N}`: A vector in which `!iszero(y2node_idx[y])`
  only if `y` is the `per` value of some node in `glo_nodes` and zero
  otherwise.
* `last_gedge_idx::E`: The highest global edge identifier already in use.
* `ppo2gbedge_idx::Array{E, 3}`: A table that translates the
  parallel-perpendicular-orientation triple to the global backward
  edge index, if it exists. If there is the need to refer to a circulation
  edge but its global identifier does not yet exists, it is created and
  saved to the table.

# Procedure

Loop through all pairs `i` and `j` (i < j) of node indexes (not
counting a possibly dummy sink node, but counting the source node),
create a edge between the two nodes if `glo_nodes[j].per - glo_nodes[i].per
<= glo_nodes[i].per`.

# Returns

1. The list of new edges, already globalized.
2. The last edge index attributed.
"""
function gen_u_fow_edges(
	glo_nodes :: Vector{Node{S, N}},
	y2node_idx :: Vector{N},
	last_gedge_idx :: E,
	ppo2gbedge_idx :: Array{E, 3}
) :: Tuple{Vector{Edge{N, E}}, E} where {S, N, E}
	u_fow_edges = Vector{Edge{N, E}}()
	isempty(glo_nodes) && return u_fow_edges, last_gedge_idx
	@assert iszero(glo_nodes[1].per)
	@assert isone(length(unique!(map(gn -> gn.par, glo_nodes))))
	@assert isone(length(unique!(map(gn -> gn.ori, glo_nodes))))
	lgei = last_gedge_idx
	for node in @view glo_nodes[2:end-1]
		bedge_ppo = (node.per, node.par, 0x03 - node.ori)
		if iszero(ppo2gbedge_idx[bedge_ppo...])
			ppo2gbedge_idx[bedge_ppo...] = (lgei += one(E))
		end
		back = ppo2gbedge_idx[bedge_ppo...]

		edge = Edge(lgei += one(E), glo_nodes[1].idx, node.idx, back)
		#(iszero(edge.head) || iszero(edge.tail)) && @show edge, @__LINE__
		push!(u_fow_edges, edge)
	end
	for i = 2:(length(glo_nodes)-1)
		for j = (i+1):length(glo_nodes)
			@assert glo_nodes[j].per > glo_nodes[i].per
			per_dist = glo_nodes[j].per - glo_nodes[i].per
			per_dist > glo_nodes[i].per && break
			iszero(y2node_idx[per_dist]) && continue
			bedge_ppo = (
				per_dist,
				glo_nodes[1].par, # all par are the same
				0x03 - glo_nodes[1].ori # all ori are the same
			)
			if iszero(ppo2gbedge_idx[bedge_ppo...])
				ppo2gbedge_idx[bedge_ppo...] = (lgei += one(E))
			end
			back = ppo2gbedge_idx[bedge_ppo...]
			edge = Edge(
				lgei += one(E), glo_nodes[i].idx, glo_nodes[j].idx, back
			)
			#(iszero(edge.head) || iszero(edge.tail)) && @show edge, @__LINE__
			push!(u_fow_edges, edge)
		end
	end
	u_fow_edges, lgei
end

"""
    gen_w_fow_edges(glo_nodes, last_gedge_idx) :: ([Edge], E)

Given a set of "globalized" nodes, and the index of the last "globalized"
edge, create the set of corresponding waste edges (already "globalized").

For now, it is not very smart and just connect the second node to the third,
third to the fourth, and so on. Ideally, it should connect every node
(except the first) to the next node *that has a backward edge associated*
(not counting the node in question, clearly, no self-loops allowed).
"""
function gen_w_fow_edges(
	glo_nodes :: Vector{Node{S, N}},
	last_gedge_idx :: E
) :: Tuple{Vector{Edge{N, E}}, E} where {S, N, E}
	w_fow_edges = Vector{Edge{N, E}}()
	lgei = last_gedge_idx
	isempty(glo_nodes) && return w_fow_edges, lgei
	@assert iszero(glo_nodes[1].per)
	@assert isone(length(unique!(map(gn -> gn.par, glo_nodes))))
	@assert isone(length(unique!(map(gn -> gn.ori, glo_nodes))))
	#sink_idx = glo_nodes[end].idx
	# Start from the second position, to avoid wasting the whole plate (if
	# the whole plate is to be wasted, this should happen at higher level).
	# Add waste arcs from every node to the next. Maybe it could be smarter
	# than this, but at least is O(n).
	for i = 2:(length(glo_nodes)-1)
		#for j = (i+1):length(glo_nodes)
		edge = Edge(lgei += 1, glo_nodes[i].idx, glo_nodes[i+1].idx, zero(E))
		#(iszero(edge.head) || iszero(edge.tail)) && @show edge, @__LINE__
		push!(w_fow_edges, edge)
		#end
	end
	w_fow_edges, lgei
end

"""
    gen_closed_flow(...)

Given some "global" structures/counters (`last_gnode_idx`, `last_gedge_idx`,
`lw2pii`, `ppo2gbedge_idx`), the piece set (`d`, `par`, `per`) and the plate
dimensions and its allowed cut orientation (`ori`, `PAR`, `PER`), this method
builds the graph that represents all allowed cuts over such plate. The only
changed parameter is `ppo2gbedge_idx`.

# Arguments

* `last_gnode_idx::N`: The highest global node identifier already in use.
* `last_gedge_idx::E`: The highest global edge identifier already in use.
* `lw2pii::AbstractArray{D, 2}`: A convenient table that translates the
  dimensions of a piece to the global piece code if there exists a piece
  that fits the description. It is not changed but a transposition of it is
  made frequently.
* `ppo2gbedge_idx::Array{E, 3}`: A table that translates the
  parallel-perpendicular-orientation triple to the global backward
  edge index, if it exists. If there is the need to refer to a circulation
  edge but its global identifier does not yet exists, it is created and
  saved to the table.
* `d::Vector{D}`: The demand of the pieces.
* `par::Vector{S}`: The size of the pieces in the dimension that is
  parallel to the cuts made in this flow.
* `per::Vector{S}`: The size of the pieces in the dimension that is
  perpendicular to the cuts made in this flow.
* `ori::UInt8`: If it is one, the flow is making vertical cuts, and
  therefore `par` is `l` and `per` is `w`. If it is two, the flow is making
  horizontal cuts and therefore `par` is `w` and `per` is `l`.
* `PAR::S`: The size of the flow/plate in the dimension that is parallel to
  the cuts made.
* `PER::S`: The size of the flow/plate in the dimension that is perpendicular
  to the cuts made.

# Returns

* A dense list of globalized nodes (all nodes in sequence).
* A sparse list of globalized nodes (if position `y` has a corresponding node,
  then `v[y] == globalized_node_id`, otherwise `v[y] == 0`).
* A dense list of the globalized edges.
* The index of the last globalized node (after the procedure).
* The index of the last globalized edge (after the procedure).
"""
function gen_closed_flow(
	last_gnode_idx :: N, last_gedge_idx :: E,
	lw2pii :: AbstractArray{D, 2}, ppo2gbedge_idx :: Array{E, 3},
	d :: Vector{D}, par :: Vector{S}, per :: Vector{S},
	ori :: UInt8, PAR :: S, PER :: S
) :: Tuple{Vector{Node{S, N}}, Vector{N}, Vector{Edge{N, E}}, N, E} where {D, S, N, E}
	# Abbreviate.
	lgni = last_gnode_idx
	lgei = last_gedge_idx

	# Reduce the pieces considered to just the ones that fit the flow/plate.
	d, par, per, n = reduce_dlw(d, par, per, PAR, PER)
	d, per = merge_duplicates(d, per)
	par = nothing # not used anymore after here
	# y2node_idx: the sparse discretization vector.
	#   If y2node_idx[y] is zero, then there is no way to reach that position y
	#   with a piece combination, otherwise, there is a way. A reachable position
	#   is first marked with one, and at the end, with the universal node index.
	# rr_fow_edges: raw restriced forward arcs. The edges (y1, y2)
	#   where y1 < y2, and both are positions inside 0:PER, and y2 - y1 is
	#   always the size of a single piece in the dimension perpendicular to
	#   the cuts.
	y2node_idx, rr_fow_edges = gen_rr_fow_edges!(
		N, d, per, PER
	)
	#@show rr_fow_edges

	dummy_sink = false
	if iszero(y2node_idx[PER])
		dummy_sink = true
		y2node_idx[PER] = one(n)
	end

	# rr_fow_edges: not raw anymore, but yet restricted
	glo_nodes, r_fow_edges, lgni, lgei = globalize!(
		y2node_idx, rr_fow_edges, lgni, lgei, lw2pii, ppo2gbedge_idx, ori, PAR
	)

	w_fow_edges, lgei = gen_w_fow_edges(glo_nodes, lgei)

	# NOTE: symmetries may be broken with extra constraints. If some edge e1 is
	# used, then no edge e2, that is both "larger" than that e1 and that may be
	# reached after traversing e1, may be used.
	u_fow_edges, lgei = gen_u_fow_edges(
		glo_nodes, y2node_idx, lgei, ppo2gbedge_idx
	)
	glo_nodes, y2node_idx, vcat(r_fow_edges, w_fow_edges, u_fow_edges), lgni, lgei
	#=
	glo_nodes, y2node_idx, append!(r_fow_edges, w_fow_edges), lgni, lgei
	=#
end

"""
    gen_nodes_and_edges(N, E, d, l, w, L, W) :: ([Node], [Edge], N, E)

Given the integer types used to number nodes and edges, the piece set, and the
original plate length, create all nodes and edges except by some dummy ones
(that are handled specially and, so, should not even be mixed here).

# Notes

* This method is aware the edge identifiers `1` to `n+4` are special.
* This method creates the special edges n+2 and n+4, but not any other in
  the `1:n+4` special identifier interval.

# High-level procedure description

* Compute the discretized lengths (widths).
* For each discretized length (width), create a flow with max width (length).
* The backward edges in not-yet-reached flow B may be referred in
  currently-being-built flow A (a subplate edge from A may enable more flow to
  a backward edge in B). To solve this problem, a structure is passed around,
  it is a matrix of three dimensions, which may indexed by the values that are
  unique for a backward edge (par_dim, cut_ori, per_dim), and stores the
  identifier for the unique backward edge with such attributes. If some forward
  edge created need to refer to a backward edge it is receives its identifier
  at that moment and is always referred by it, the corresponging `Edge` struct
  however, is only created at the final of the procedure, when all backward
  edges are known.
* Backward Edges are similar to waste edges in the fact both of they have the
  `back` field as zero (they never enable another backward edge), however only
  backward edges have tails pointing to the source nodes of flows (and,
  consequently, the only ones with tails before heads).
"""
function gen_all_edges(
	::Type{N}, ::Type{E},
	d :: Vector{D}, l :: Vector{S}, w :: Vector{S},
	L :: S, W :: S
) :: Tuple{Vector{Node{S, N}}, Vector{Edge{N, E}}, N, E} where {D, S, N, E}
	@assert length(d) == length(l)
	@assert length(d) == length(w)
	n = length(d)
	lgni = zero(N)
	lgei = convert(E, n + 4) # the first n + 4 edge indexes are reserved
	edges = Vector{Edge{N, E}}()
	nodes = Vector{Node{S, N}}()
	lw2pii = fill(zero(E), L, W)
	for pii = 1:n
		(l[pii] > L || w[pii] > W) && continue
		if !iszero(lw2pii[l[pii], w[pii]])
			@warn(
				"CAUTION: the instance has two pieces with the exact same length" *
				" and width, the first is piece nº $(lw2pii[l[pii], w[pii]]) and" *
				" the second is nº $(pii), both have length $(l[pii]) and width " *
				"$(w[pii])"
			)
		end
		lw2pii[l[pii], w[pii]] = pii
	end
	max_LW = max(L, W)
	ppo2gbedge_idx = fill(zero(E), max_LW, max_LW, 2)
	origin_vf_by_L = fill(zero(N), L)
	origin_hf_by_W = fill(zero(N), W)
	vflows_by_L = [Vector{N}() for _ = 1:L]
	hflows_by_W = [Vector{N}() for _ = 1:W]

	# Mark the two root backward edges to be created at the end
	ppo2gbedge_idx[L, W, 0x01] = n + 2
	ppo2gbedge_idx[W, L, 0x02] = n + 4

	# Create the cutting of the root plate with the default orientation
	# (vertical).
	vroot_nodes, vflows_by_L[L], vroot_edges, lgni, lgei = gen_closed_flow(
		lgni, lgei, lw2pii, ppo2gbedge_idx, d, l, w, 0x01, L, W
	)
	@assert !isempty(vroot_nodes)
	@assert iszero(vroot_nodes[1].per)
	@assert vroot_nodes[1].par == L
	origin_vf_by_L[L] = vroot_nodes[1].idx
	disc_W = map(n -> n.per, vroot_nodes)
	#@show vroot_edges
	append!(edges, vroot_edges)
	append!(nodes, vroot_nodes)

	#for ppo in CartesianIndices(ppo2gbedge_idx) # DEBUG PRINT FOR
	#  iszero(ppo2gbedge_idx[ppo]) && continue
	#  @show ppo
	#  @show ppo2gbedge_idx[ppo]
	#end

	# Create the forward edge that inverts the default orientation
	# of the root plate.
	push!(edges, Edge{N, E}(
		n + 3, origin_vf_by_L[L], lgni, n + 4
	))
	#@show last(edges), @__LINE__

	hroot_nodes, hflows_by_W[W], hroot_edges, lgni, lgei = gen_closed_flow(
		lgni, lgei, lw2pii', ppo2gbedge_idx, d, w, l, 0x02, W, L
	)
	origin_hf_by_W[W] = hroot_nodes[1].idx
	disc_L = map(n -> n.per, hroot_nodes)
	#@show hroot_edges
	append!(edges, hroot_edges)
	append!(nodes, hroot_nodes)

	#for ppo in CartesianIndices(ppo2gbedge_idx) # DEBUG PRINT FOR
	#  iszero(ppo2gbedge_idx[ppo]) && continue
	#  @show ppo
	#  @show ppo2gbedge_idx[ppo]
	#end

	#unique_lw
	#@show disc_L
	#@show disc_W
	for L_ in @view disc_L[2:end-1]#setdiff(unique(l), [L])#
		vnodes, vflows_by_L[L_], vedges, lgni, lgei = gen_closed_flow(
			lgni, lgei, lw2pii, ppo2gbedge_idx, d, l, w, 0x01, L_, W
		)
		origin_vf_by_L[L_] = vnodes[1].idx
		append!(edges, vedges)
		append!(nodes, vnodes)
	end
	for W_ in @view disc_W[2:end-1]#setdiff(unique(w), [W])#
		hnodes, hflows_by_W[W_], hedges, lgni, lgei = gen_closed_flow(
			lgni, lgei, lw2pii', ppo2gbedge_idx, d, w, l, 0x02, W_, L
		)
		origin_hf_by_W[W_] = hnodes[1].idx
		append!(edges, hedges)
		append!(nodes, hnodes)
	end
	#@show vflows_by_L
	#@show hflows_by_W

	for ppo in CartesianIndices(ppo2gbedge_idx)
		pari, peri, orii = ppo[1], ppo[2], ppo[3]
		iszero(ppo2gbedge_idx[ppo]) && continue
		#@show ppo
		if orii == 0x01
			while !iszero(peri) && iszero(vflows_by_L[pari][peri]); peri -= 1; end
			@assert !iszero(peri)
			push!(edges, Edge{N, E}(
				ppo2gbedge_idx[ppo], vflows_by_L[pari][peri], origin_vf_by_L[pari],
				zero(N)
			))
			edge = last(edges)
			#(iszero(edge.head) || iszero(edge.tail)) && @show edge, ppo, peri, @__LINE__
		else
			@assert orii == 0x02
			while !iszero(peri) && iszero(hflows_by_W[pari][peri]); peri -= 1; end
			push!(edges, Edge{N, E}(
				ppo2gbedge_idx[ppo], hflows_by_W[pari][peri], origin_hf_by_W[pari],
				zero(N)
			))
			edge = last(edges)
			#(iszero(edge.head) || iszero(edge.tail)) && @show edge, ppo, peri, @__LINE__
		end
	end

	nodes, edges, lgni, lgei
end

end # module

