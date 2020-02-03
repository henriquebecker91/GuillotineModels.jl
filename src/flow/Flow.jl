module Flow

module Args
	import ...Utilities
	import ...Utilities.Args: Arg
	#export accepted_arg_list, throw_if_incompatible_options
	# No extra flags for this model yet.
	Utilities.Args.accepted_arg_list(::Val{:Flow}) = Vector{Arg}()
	function Utilities.Args.throw_if_incompatible_options(
		::Val{:Flow}, p_args
	)
		nothing
	end
end

# Used for building the model, used by module Format below.
include("Enumeration.jl")
using .Enumeration

# Just include it to make it available to users of this module.
include("Format.jl")

using JuMP

import ..build_model
"""
    build_model(::Val{:Flow}, model, d, p, l, w, L, W, [options])

The Flow model is basically the PPG2KP model but using flow constraints.
The result is a considerably larger model with a slight enhancement
on the relaxation bound.

# Arguments

The required arguments are the same as all `build_model` methods, and
the optional argument `options` is always ignored.

# High-Level Model Description

Both length and width of the original plate may be seen as line segments.
From these line segments, let us just consider a subset of the points:
the first point (zero), the last point (L or W), and the ones in the middle
that match a normal cut position.

For each width cut position, there is an instance of the length point set
representing plates with all possible lengths but a fixed width.
The analogue (each width cut position ... length point set) also exists.
Each of these point set instances is the vertex set of a strongly connected
directed graph.
There are no edges connecting one of such graphs with another (they interact by
special constrainsts that cannot be seen as 'flow constraints'/edges).

The variables of the model are the edges of those graphs.
The edges are may be divided into forward edges and backward/circulation
edges.
The forward edges go from a node representing a point closer to the start
of the line segment to a node representing a point farther from the start.
The forward edges are all edges in which the origin is closer to zero
(i.e., the first point of that set) than the destination point.
The backward/circulation edges always go from a non-first/zero point to
the first/zero point.
The forward edges represent cutting part of the plate for use or trimming.
The backard edges represent the availibity of a plate of some dimensions
to be cut in some orientation.
For example, if there is two units of flow flowing by a forward edge that
goes from point length 15 to point length 35 of the set with fixed width
50, then there are at least two plates of width 50 that have a subplate
of length 20 (and width 50) cut from two parallel cuts, one in position
15 and another in 35. The two parent plates may have the same dimensions
(i.e., their length is also equal to each other), or they may have
different total length (they only share the same width); for an example,
one of them may stop at length 35 and the other may have length 1000.
Considering this same example, if there are a plate of size 35x50 and
another of size 1000x50 available, then there is one of flow in the
backward edge from point length 35 to point zero and from point length
1000 to point zero (all inside the set with fixed width 50).

Finally, forward edges may be subdivided into waste edges, piece edges and
subplate edges.
Waste edges do nothing besides letting the flow pass through them.
The variables that represent the flow amount of a piece edge are present in
the objective function (being multiplied by the profit of the corresponding
piece).
Subplate edges "introduce" new flow in the system, they do not only allow the
flow to pass through them but also increase the flow in a correspondent
circulation edge of the (other) graph that represents the plate they are
extracting, but to be cut in the opposite orientation.
For example, an increase in the flow of the subplate edge from length 50 to
length 150 of the set with fixed width 200 will increase by the same amount the
flow of the backward edge from width 200 to width zero of the set with fixed
length 100.

The number and dimension of the original plates available is given by all
backward edges that have their value set as non-zero a priori.

# Slightly Outdated and Misleading Low-Level Explanation

## Variables

* `cut_edges[PP, PO, CP]`: Integer. The cutting edges. An arc pertains to a
   flow PP (Parent Plate), its start point is PO (POint or POsition), and it
   allows one more flow unit pass through CP (Child Plate). If the child plate
   is zero, then the arc does not have a child plate. Arcs that do not have a
   child are either: (1) waste edges that connect PO to the end of PP; (2)
   piece/"final plate" edges that have a nonzero coefficient in the value of
   the objective function; (3) backward edges that return the flow from sink
   to source.
* `cir_edges[P]`: Integer. The circulation edges. A circulation arc pertains
   to some plate P, and connect the end of the plate/flow back to its
   beggining.

## Objective function:

Maximize the profit of the pieces cut.

`sum(p[CP] * cut_edges[_, _, CP is a final plate])`

## Constraints

* There is exactly one of the original plate, so its circulation arc is
  available since start, and cannot ever go over one. It must be one.
  `sum(cir_edges[P is the original plate]) == 1`
* The flow of a plate must go round.
  - `sum(cut_edges[PP, PO1, CP]) == sum(cut_edges[PP, PO2, vary CP]) :
    for all PP, PO1 + CP length == PO2`
  - `cir_edges[PP] == sum(cut_edges[PP, PO is the end of PP, vary CP]) :
      for all PP`
  - `cir_edges[PP] == sum(cut_edges[PP, PO is the start of PP, vary CP]) :
      for all PP`
* If you cut a subplate anywhere then you allow its flow to circle again.
  - `sum(cir_edges[P]) <= sum(cut_edges[vary PP, vary PO, P]) :
    for all P that are intermediary plates`
* The number of pieces of some type is always less than or equal to the demand.
  - `sum(cut_edges[vary PP, vary PO, CP]) <= d[CP] :
      for all CP that are final plates (i.e., pieces)`

## Conventions:
* 0 - The plate is waste. This may only appear at the CP index of cut_edges.
* 1:n - The plate is a final plate corresponding to a piece index.
* n+1 - A dummy edge that is upper bounded at one and allows one extra flow
  to the original plate default cutting (vertical) backward edge (i.e.,
  n+2). Necessary to avoid n+2 to be restricted by a possible edge that
  cuts the whole original plate alternative cutting (horizontal). In other
  words, without this, if a piece has the same length as the plate, then
  the model has a deadlock problem, and will give a zero objective problem.
* n+2 - The original plate default cutting (vertical) backward edge.
* n+3 - The original plate vertical forward edge that takes the whole plate
  and allows it to be horizontally (i.e., allow more flow to n+3).
* n+4 - The original plate alternative cutting (horizontal) backward edge.
  It is allowed by cutting n+2.
* n+5:m - The remaining intermediary plates.
"""
function build_model(
	::Val{:Flow}, model, d :: Vector{D}, p :: Vector{P},
	l :: Vector{S}, w :: Vector{S}, L :: S, W :: S,
	options :: Dict{String, Any} = Dict{String, Any}()
) :: Tuple{Vector{Node{S, P}}, Vector{Edge{P, P}}} where {D, S, P}
	# TODO: fix (L, W, l, w, d) = (20, 5, [20, 4, 18], [4, 5, 2], [1, 1, 1])
	@assert length(d) == length(l)
	@assert length(d) == length(w)
	@assert length(d) == length(p)
	num_piece_types = convert(D, length(d))

	N = E = P
	# TODO: check if N and E will be bubbled up.
	# NOTE: the profits are not needed, the first num_piece_types edges will
	# be dummy edges that have the profit value in the objective function.
	nodes_data, edges_data, last_gnode_idx, last_gedge_idx = gen_all_edges(
		N, E, d, l, w, L, W
	)
	#@show edges_data

	head2indxs = [Vector{E}() for _ = 1:last_gnode_idx]
	tail2indxs = [Vector{E}() for _ = 1:last_gnode_idx]
	back2indxs = [Vector{E}() for _ = 1:last_gedge_idx]
	for e in edges_data
		@assert e.indx <= last_gedge_idx
		@assert e.head <= last_gnode_idx
		iszero(e.head) && @show e
		@assert e.tail <= last_gnode_idx
		iszero(e.tail) && @show e
		push!(head2indxs[e.head], e.indx)
		push!(tail2indxs[e.tail], e.indx)
		!iszero(e.back) && push!(back2indxs[e.back], e.indx)
	end
	# Create the dummy edge that enables the "original plate default orientation
	# backward edge" to pass one unit of flow.
	dummy_vroot_enabler = num_piece_types + one(D)
	vroot_back_edge = dummy_vroot_enabler + one(D)
	#@assert isempty(back2indxs[vroot_back_edge])
	push!(back2indxs[vroot_back_edge], dummy_vroot_enabler)

	@variable(model, edge[1:last_gedge_idx] >= 0, Int)

	@objective(model, Max,
		sum(p[pii] * edge[pii] for pii = 1:num_piece_types)
	)

	# There is just one of the original plate. Sets the dummy edge to one and
	# consequently allow one unit of flow to the original plate default
	# orientation
	# backward edge.
	@constraint(model, edge[dummy_vroot_enabler] == 1)

	for i = 1:last_gnode_idx
		if isempty(head2indxs[i]) || isempty(tail2indxs[i])
			@show i
			@show head2indxs[i]
			@show tail2indxs[i]
		end
		#@assert !isempty(head2indxs[i]) && !isempty(tail2indxs[i])
		@constraint(model,
			sum(edge[head2indxs[i]]) == sum(edge[tail2indxs[i]])
		)
	end

	for i = 1:num_piece_types # respect demand
		@constraint(model, edge[i] <= d[i])
	end

	# Not all edges are cir_edges, however, we do not have a separate list at
	# the moment. The easiest way is to traverse all edges, if an edge has
	# a non-empty back2indxs, then it has other edges allowing more flow to it
	# and, consequently, it is a backward/circulation edge.
	for i = 1:last_gedge_idx
		@assert i ∉ back2indxs[i]
		if !isempty(back2indxs[i])
			@constraint(model,
				edge[i] <= sum(edge[back2indxs[i]])
			)
		end
	end

	# Create fake edges for the pieces.
	for i = 1:num_piece_types
		push!(edges_data, Edge{N, E}(i, zero(N), zero(N), zero(E)))
	end
	# Create the dummy edge for the number of original plates.
	push!(edges_data, Edge{N, E}(
		num_piece_types + 1, zero(N), zero(N), num_piece_types + 2
	))
	# Order the edges by their index.
	sort!(edges_data, by = e -> e.indx)

	nodes_data, edges_data
end # build_model_no_symmbreak

end # module

