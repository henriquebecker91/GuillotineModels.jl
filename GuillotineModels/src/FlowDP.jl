module FlowDP

export Node, Edge, FlowModel, gen_all_edges

struct Node{S, N}
  idx :: N # The node index in the universe of nodes.
  par :: S # 
  per :: S # 
  ori :: UInt8 # 
end

struct Edge{N, E}
  indx :: E # The indx index in the universe of edges.
  head :: N # The node from where the flow comes. Source.
  tail :: N # The node from where the flow ends. Sink.
  back :: E # The backward edge global index, if it exists, otherwise zero. 
end

#=
struct FlowModel{S, N, E}
  nodes :: Vector{Node{S, N}}
  edges :: Vector{Edge{N, E}}
end

function append!(
  lhs :: FlowModel{S, N, E}, rhs :: FlowModel{S, N, E}
) :: FlowModel{S, N, E}
  append!(lhs.nodes, rhs.nodes)
  append!(lhs.edges, rhs.edges)
end
=#

#     gen_rr_fow_edges!
# Take the pieces demand, the pieces length (or width), and the plate
# corresponding dimension. Returns: (1) a vector with the same size as the
# plate dimension where vector[y] == 1 if there is a linear combination
# of piece sizes that give `y` and zero otherwise; (2) every pair of such
# "one-positions" for which there exists a piece with the same exact size as
# the distance between the two positions.
# NOTE: N is the type of the node indexes, to allow the method to already
# create one of the returned structures with the right element type.
# NOTE: THIS ASSUMES ALL PIECES FIT THE PLATE ON BOTH DIMENSIONS.
function gen_rr_fow_edges!(
  ::Type{N}, # eltype of the marked vector
  d :: Vector{D}, l :: Vector{S}, L :: S
) :: Tuple{Vector{N}, Vector{Tuple{S, S}}} where {D, S, N}
  @assert length(d) == length(l)
  n = length(d)
  marked :: Vector{N} = fill(zero(N), L)
  edges = Vector{Tuple{S, S}}() # start point and end point
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
      push!(edges, (y_, y += l[1]))
      marked[y] = one(N)
    end
  end

  marked, edges
end

# TODO: decide where the pieces will be pieces will be reduced, and if
# it is necessary to expand them back, and how this will be done. The
# first version of the code will not reduce them and instead will just 
# check the dimensions every time.
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

# TODO: this could be O(n log n) instead of pseudopolynomial
function merge_duplicates(
  d :: Vector{D}, per :: Vector{S}
) :: Tuple{Vector{D}, Vector{S}} where {D, S}
  @assert length(d) == length(per)
  #isempty(d) && return d, per

  y2d = fill(zero(D), maximum(per))
  for i in eachindex(d)
    y2d[per[i]] += d[i] 
  end

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

# Takes the y2node_idx discretization and globalize it (i.e., change the ones
# to the global node indexes). Create a list of discretized positions as sub-
# product. Create a list of globalized forward edges based on the list of
# raw_fow_edges passed as argument. Also return the updated values of
# last_gnode_idx and last_gedge_idx after registering the edges and nodes.
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
    (iszero(head) || iszero(tail)) && @show edge, @__LINE__
    push!(glo_fow_edges, edge)
    # If this foward edge cuts a plate with the exact same size as some piece
    # then we add the "final plate"/piece edge too.
    yd = y2 - y1
    max_par = size(lw2pii, 1)
    max_per = size(lw2pii, 2)
    if PAR <= max_par && yd <= max_per && !iszero(lw2pii[PAR, yd])
      indx = (lgei += one(E))
      edge = Edge{N, E}(indx, head, tail, lw2pii[PAR, yd])
      push!(glo_fow_edges, edge)
    end
  end

  glo_nodes, glo_fow_edges, lgni, lgei
end

# TODO: this is generating all unrestricted arcs again, solve this
function gen_u_fow_edges(
  glo_nodes :: Vector{Node{S, N}},
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

    edge = Edge(lgei += 1, glo_nodes[1].idx, node.idx, back)
    (iszero(edge.head) || iszero(edge.tail)) && @show edge, @__LINE__
    push!(u_fow_edges, edge)
  end
  for i = 2:(length(glo_nodes)-1)
    for j = (i+1):length(glo_nodes)
      @assert glo_nodes[j].per > glo_nodes[i].per
      bedge_ppo = (
        glo_nodes[j].per - glo_nodes[i].per,
        glo_nodes[1].par, # all par are the same
        0x03 - glo_nodes[1].ori # all ori are the same
      )
      if iszero(ppo2gbedge_idx[bedge_ppo...])
        ppo2gbedge_idx[bedge_ppo...] = (lgei += one(E))
      end
      back = ppo2gbedge_idx[bedge_ppo...]

      edge = Edge(
        # TODO: this does not work, those are waste edges, not use edges
        lgei += 1, glo_nodes[i].idx, glo_nodes[j].idx, back
      )
      (iszero(edge.head) || iszero(edge.tail)) && @show edge, @__LINE__
      push!(u_fow_edges, edge)
    end
  end
  u_fow_edges, lgei
end

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
  sink_idx = glo_nodes[end].idx
  for node = @view glo_nodes[2:end-1]
    edge = Edge(lgei += 1, node.idx, sink_idx, zero(E))
    (iszero(edge.head) || iszero(edge.tail)) && @show edge, @__LINE__
    push!(w_fow_edges, edge)
  end
  w_fow_edges, lgei
end

# TODO: We need a better way to deal with the pieces. Technically, two
#   distinct pieces may have the same length and width, but cannot be
#   merged because they have different profits. For now, we will ignore
#   this case, as this case does not seem common (may be even prohibited
#   in some instance generations), but in the future, we cannot identify
#   single pieces by their dimensions, we need to repeat whatever is being
#   done for each piece with the same dimensions, or pass a true identifier
#   along.
function gen_closed_flow(
  last_gnode_idx :: N, # The highest global node identifier already in use.
  last_gedge_idx :: E, # The highest global edge identifier already in use.
  lw2pii :: AbstractArray{D, 2}, # A convenient table that translates the
    # dimensions of a piece to the global piece code if there exists a piece
    # that fits the description. It is not changed but a transposition of it is
    # made frequently.
  ppo2gbedge_idx :: Array{E, 3}, # A table that translates the
    # parallel-perpendicular-orientation triple to the global backward
    # edge index, if it exists. If there is the need to refer to a circulation
    # edge but its global identifier does not yet exists, it is created and
    # saved to the table.
  d :: Vector{D}, # The demand of the pieces.
  par :: Vector{S}, # The size of the pieces in the dimension that is parallel
    # to the cuts made in this flow.
  per :: Vector{S}, # The size of the pieces in the dimension that is
    # perpendicular to the cuts made in this flow.
  ori :: UInt8, # If it is one, the flow is making vertical cuts, and therefore
    # `par` is `l` and `per` is `w`. If it is two, the flow is making
    # horizontal cuts and therefore `par` is `w` and `per` is `l`.
  PAR :: S, # The size of the flow/plate in the dimension that is parallel to
    # the cuts made.
  PER :: S # The size of the flow/plate in the dimension that is perpendicular
    # to the cuts made.
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

  # TODO: make the first version without the unrestricted edges and check if
  # it get the same values as the old_bender. Execute against A5 and
  # check if it get the restricted value, then implement the unrestricted
  # edges and check if it get the unrestricted value for A5.
  # The unrestrict edges are generated by the following procedure:
  # loop through all pairs `i` and `j` (i < j) of disc_per indexes (not
  # counting a possibly dummy sink node, but counting the source node),
  # create a edge between the two nodes if `disc_per[j] - disc_per[i] <=
  # disc_per[i]`.
  # NOTE: symmetries may be broken with extra constraints. If some edge e1 is
  # used, then no edge e2, that is both "larger" than that e1 and that may be
  # reached after traversing e1, may be used.
  #=
  u_fow_edges, lgei = gen_u_fow_edges(
    glo_nodes, lgei, ppo2gbedge_idx
  )
  glo_nodes, y2node_idx, vcat(r_fow_edges, w_fow_edges, u_fow_edges), lgni, lgei
  =#
  glo_nodes, y2node_idx, append!(r_fow_edges, w_fow_edges), lgni, lgei
end

#= Completely inutile. We need to be able to traverse the per dimension.
function add_nodes_to_dict(
  d :: Dict{Tuple{S, S, UInt8}, N},
  nodes :: Vector{Node{S, N}}
) :: Dict{Tuple{S, S, UInt8}, N}
  for node in nodes
    ppo = (node.par, node.per, node.ori)
    @assert !haskey(d, ppo)
    d[ppo] = node.idx
  end

  d
end
=#

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
    @assert iszero(lw2pii[l[pii], w[pii]])
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
  
  for L_ in @view disc_L[2:end-1]#setdiff(unique(l), L)#
    vnodes, vflows_by_L[L_], vedges, lgni, lgei = gen_closed_flow(
      lgni, lgei, lw2pii, ppo2gbedge_idx, d, l, w, 0x01, L_, W
    )
    origin_vf_by_L[L_] = vnodes[1].idx
    append!(edges, vedges)
    append!(nodes, vnodes)
  end
  for W_ in @view disc_W[2:end-1]#setdiff(unique(w), W)#
    hnodes, hflows_by_W[W_], hedges, lgni, lgei = gen_closed_flow(
      lgni, lgei, lw2pii', ppo2gbedge_idx, d, w, l, 0x02, W_, L
    )
    origin_hf_by_W[W_] = hnodes[1].idx
    append!(edges, hedges)
    append!(nodes, hnodes)
  end

  # TO CREATE THE CIRCULATION EDGES AT THE END, WE NEED TO KNOW THE NAME OF
  # THE NODES FOR THE ORIGIN OF EVERY FLOW WITH PAR WITH SOME SIZE, AND THE
  # INTERNAL NODE 
  #@show vflows_by_L
  #@show origin_vf_by_L
  #@show hflows_by_W
  #@show origin_hf_by_W
  for ppo in CartesianIndices(ppo2gbedge_idx)
    pari, peri, orii = ppo[1], ppo[2], ppo[3]
    iszero(ppo2gbedge_idx[ppo]) && continue
    @show ppo
    if orii == 0x01
      while !iszero(peri) && iszero(vflows_by_L[pari][peri]); peri -= 1; end
      @assert !iszero(peri)
      push!(edges, Edge{N, E}(
        ppo2gbedge_idx[ppo], vflows_by_L[pari][peri], origin_vf_by_L[pari],
        zero(N)
      ))
      edge = last(edges)
      (iszero(edge.head) || iszero(edge.tail)) && @show edge, ppo, peri, @__LINE__
    else
      @assert orii == 0x02
      while !iszero(peri) && iszero(hflows_by_W[pari][peri]); peri -= 1; end
      push!(edges, Edge{N, E}(
        ppo2gbedge_idx[ppo], hflows_by_W[pari][peri], origin_hf_by_W[pari],
        zero(N)
      ))
      edge = last(edges)
      (iszero(edge.head) || iszero(edge.tail)) && @show edge, ppo, peri, @__LINE__
    end
  end
  
  # There is a flow for every discretized point on the dimension parallel to
  # the cut, for each of the dimensions/'cut orientations'. Therefore, the
  # worst case is L + W.
  #
  # PSEUDOCODE (first version, not very optimized)
  # compute the discretized lengths (widths)
  # for each discretized length (width), create a flow with max width (length)
  # the circulation arcs connect a position in some flow to the source of the
  #   same flow, they may be identified by which dimension is parallel to
  #   the cut, the size of that dimension (these two identify the flow), and
  #   the size of the perpendicular dimension (this one identifies which one
  #   of the multiple circulation arcs inside a flow is being referred to)
  # the circulation arcs are referred between different closed flows, so their
  #   global name needs to be created at the moment some arc has the
  #   circulation arc as their CP field
  # a global structure is passed around, it is a matrix of three dimension,
  #   which may indexed by the values that are unique for a circulation arc
  #   (par_dim, cut_ori, per_dim), and stores the code for the unique
  #   circulation arc with such attributes. All the circulation arc variables
  #   are created at the end. They will have CP zero as would a waste arc,
  #   because a circulation arc never enables another circulation arc, and
  #   they are the only arcs for which the first vertex comes after the second
  #   (and the second is always the source), but none of those peculiarities
  #   prevent us of modeling the same way waste arcs are modeled.
  #
  nodes, edges, lgni, lgei
end

end # module

