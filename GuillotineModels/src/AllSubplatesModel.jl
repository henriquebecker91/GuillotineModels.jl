module AllSubplatesModel

using JuMP
include("GuillotinePlatesDP.jl")
using .GuillotinePlatesDP

# HIGH LEVEL EXPLANATION OF THE MODEL
#
# Variables:
#
# `picuts[n, pii]`: Integer. The number of pieces `pii` generated from
#   subplates of type `n`.
# `cuts_made[n1, n2, n3]`: Integer. The number of subplates of type
#   `n1` that are cut into subplates `n2` and `n3` (horizontal and
#   vertical cuts are together for now).
#
# Objective function:
#
# Maximize the profit of the pieces cut.
#   sum(p[pii] * picuts[_, pii])
#
# Constraints:
#
# There is exactly one of the original plate, which may be used for cutting
# or extracting a piece.
#   sum(picuts[1, _]) + sum(cuts_made[1, _,  _]) <= 1
# The number of subplates available depends on the number of plates that have
# it as children.
#   sum(picuts[n1>1, _]) + sum(cuts_made[n1>1, _, _]) <=
#     sum(cuts_made[_, n2, n3])
#     where n2 == n1 or n3 == n1, doubling cuts_made[_, n2, n3] if n2 == n3
# The number of pieces of some type is always less than or equal to the demand.
#   sum(picuts[_, pii]) <= d[pii]
#
# Unnecessary constraints:
#
# The number of times a pair pli-pii appear is at most the min between: 1)
# d[pii] and the number of subplates pli that fit in the original plate.
#   sum(picuts[n, pii]) <= min(d[pii], max_fits[n])
function build_model_no_symmbreak(
  model, d :: Vector{D}, p :: Vector{P}, l :: Vector{S}, w :: Vector{S},
  L :: S, W :: S; only_binary = false
) where {D, S, P}
  @assert length(d) == length(l) && length(l) == length(w)
  num_piece_types = convert(D, length(d))

  sllw = SortedLinkedLW(D, l, w)
  pli_lwb, hcuts, vcuts, np = gen_cuts(P, d, sllw, L, W)
  num_plate_types = length(pli_lwb)
  hvcuts = vcat(hcuts, vcuts)

  # Order the plate-piece pairs by my ranking of importance: how much
  # absolute are is wasted. Other rankings include: how much relative area
  # is wasted; how much is the profit of the a squared unit of the plate used;
  # group them by pii just for making it easier to the demand constraint.
  #=sort!(np, lt = function((pli1, pii1), (pli2, pii2))
    pli1l, pli1w, _ = pli_lwb[pli1]
    pli2l, pli2w, _ = pli_lwb[pli2]
    pii1l, pii1w = l[pii1], w[pii1]
    pii2l, pii2w = l[pii2], w[pii2]
    convert(P, pli1l - pii1l) * (pli1w - pii1w) < (
    convert(P, pli2l - pii2l) * (pli2w - pii2w))
  end)=#

  # pli2pair: inverse index which given a plate index will return a list of
  # all picuts indexes (np variable) that cut some piece from some plate.
  pli2pair = [Vector{P}() for _ = 1:num_plate_types]
  # pii2pair: the same as pli2pair but for piece indexes.
  pii2pair = [Vector{P}() for _ = 1:num_piece_types]

  # The vectors below all have the same length as the number of plate types (to
  # allow indexing by plate type). The value of a position is a vector of
  # arbitrary length and irrelevant index, the values of this inner vector are
  # cut indexes. Such cut indexes are related to the plate that is the index of
  # the outer vector.
  # any CHILD plate to respective CUT indexes
  child2cut = [Vector{P}() for _ = 1:num_plate_types]
  # any PARENT plate to respective CUT indexes
  parent2cut = [Vector{P}() for _ = 1:num_plate_types]

  # Initialize all inverse indexes.
  for i in eachindex(hvcuts)
    parent, fchild, schild = hvcuts[i]
    push!(parent2cut[parent], i)
    push!(child2cut[fchild], i)
    @assert !iszero(schild)
    push!(child2cut[schild], i)
  end

  for i in eachindex(np)
    pli, pii = np[i]
    push!(pli2pair[pli], i)
    push!(pii2pair[pii], i)
  end

  # If all pieces have demand one, a binary variable will suffice to make the
  # connection between a piece type and the plate it is extracted from.
  naturally_only_binary = all(di -> di <= 1, d)
  if naturally_only_binary || only_binary
    # only_binary is equal to naturally_only_binary for now, but the idea is
    # that only_binary will expand the number of binary variables to account
    # for picuts that can repeat (plates that can appear more than one time
    # and that may have the same piece extracted from them)
    @variable(model, picuts[1:length(np)], Bin)
  else
    @variable(model, picuts[1:length(np)] >= 0, Int)
    #@variable(model,
    #  0 <= picuts[i = 1:length(np)] <=
    #    min(pli_lwb[np[i][1]][3], d[np[i][2]]),
    #Int)
  end

  if only_binary
    @variable(model, cuts_made[1:length(hvcuts)], Bin)
  else
    @variable(model, cuts_made[1:length(hvcuts)] >= 0, Int)
  end

  # The objective function is to maximize the profit made by extracting
  # pieces from subplates.
  @objective(model, Max,
    sum(p[pii] * sum(picuts[pii2pair[pii]]) for pii = 1:num_piece_types)
  )

  # ub0: trivial area upper bound
  #@constraint(model, sum(a[j] * x[i, j] for i = 1:N, j = 1:T) <= L*W)

  # c1: There is just one of the original plate, and so it can be only used
  # to extract a single piece xor make a single cut that would make two new
  # subplates available.
  @constraint(model,
    sum(picuts[pli2pair[1]]) + sum(cuts_made[parent2cut[1]]) <= 1
  )

  # c2: for each subplate type that is not the original plate, such subplate
  # type will be available the number of times it was the child of a cut,
  # subtracted the number of times it had a piece extracted or used for
  # further cutting.
  for pli in 2:num_plate_types
    @constraint(model,
      sum(picuts[pli2pair[pli]]) + sum(cuts_made[parent2cut[pli]]) <=
      sum(cuts_made[child2cut[pli]])
    )
  end

  # c2.5: The amount of each subplate type generated by cuts is bounded by the
  # amount that can be cut from the original plate.
  for pli in 2:num_plate_types
    @constraint(model,
      sum(cuts_made[parent2cut[pli]]) <= pli_lwb[pli][3]
    )
  end

  # c3: the amount of each piece type extracted from different plate types
  # cannot surpass the demand for that piece type.
  for pii in 1:num_piece_types
    @constraint(model, sum(picuts[pii2pair[pii]]) <= d[pii])
  end

  model, hvcuts, pli_lwb, np
end # build_model_no_symmbreak

# HIGH LEVEL EXPLANATION OF THE MODEL
#
# Variables:
#
# `pieces_sold[pii]`: Integer. The number of pieces `pii` sold. Is the minimum
#   between the demand of the piece and the amount of plates generated and not
#   used that have exactly the same size as the piece.
# `cuts_made[n1, n2, n3]`: Integer. The number of subplates of type `n1` that
#   are cut into subplates `n2` and `n3` (horizontal and vertical cuts are
#   together for now). As this is a symmetry-breaking model, the plate types
#   are not only each distinct `l` and `w` but each different `l`, `w`, and
#   `symm` (that marks if the plate can be cut only horizontally, only
#   vertically, or both ways).
#
# Objective function:
#
# Maximize the profit of the pieces sold.
#   sum(p[pii] * pieces_sold[pii])
#
# Constraints:
#
# There is exactly one of the original plate, which may be used for cutting
# or extracting a piece.
#   sum(pieces_sold[plates exactly the size of the original plate]) +
#   sum(cuts_made[1, _,  _]) <= 1
# The number of subplates available depends on the number of plates that have
# it as children.
#   sum(cuts_made[n1>1, _, _]) <= sum(cuts_made[_, n2, n3])
#     where n2 == n1 or n3 == n1, doubling cuts_made[_, n2, n3] if n2 == n3
# The number of pieces sold is bounded both by the demand of the piece type and
# the the number of unused plates with the same size as the piece.
#   sum(pieces_sold[pii]) <= d[pii]
#   sum(pieces_sold[pii]) <= cuts_made[_, n2, n3] - cuts_made[n2 or n3, _, _]
#     where n2 or n3 has the same size as pii, fixing for when n2 == n3
#
# Unnecessary constraints:
#
# The amount of times a plate type may be cut is bounded by how many of them
# could fit the original plate. Note that we ignore the symmetry tag here
# and group all the plates with the same `l` and `w` but distinct symmetry tag.
#   sum(cuts_made[plates sharing `l` and `w`, _, _]) <= (L ÷ l) * (W ÷ w)
function build_model_with_symmbreak(
  model, d :: Vector{D}, p :: Vector{P}, l :: Vector{S}, w :: Vector{S},
  L :: S, W :: S; only_binary = false
) where {D, S, P}
  @assert length(d) == length(l) && length(l) == length(w)
  num_piece_types = convert(D, length(d))

  sllw = SortedLinkedLW(D, l, w)
  pli2lwsb, hcuts, vcuts, pii2plis, pli2piis, same_size_plis =
    gen_cuts_sb(P, d, sllw, L, W)
  num_plate_types = length(pli2lwsb)
  hvcuts = vcat(hcuts, vcuts)

  # The vectors below all have the same length as the number of plate types (to
  # allow indexing by plate type). The value of a position is a vector of
  # arbitrary length and irrelevant index, the values of this inner vector are
  # cut indexes. Such cut indexes are related to the plate that is the index of
  # the outer vector.
  # any CHILD plate to respective CUT indexes
  child2cut = [Vector{P}() for _ = 1:num_plate_types]
  # any PARENT plate to respective CUT indexes
  parent2cut = [Vector{P}() for _ = 1:num_plate_types]

  # Initialize all inverse indexes.
  for i in eachindex(hvcuts)
    parent, fchild, schild = hvcuts[i]
    push!(parent2cut[parent], i)
    push!(child2cut[fchild], i)
    !iszero(schild) && push!(child2cut[schild], i)
  end

  # If all pieces have demand one, a binary variable will suffice to make the
  # connection between a piece type and the plate it is extracted from.
  naturally_only_binary = all(di -> di <= 1, d)
  if naturally_only_binary || only_binary
    # only_binary is equal to naturally_only_binary for now, but the idea is
    # that only_binary will expand the number of binary variables to account
    # for piece solds that can repeat (i.e., have demand more than one).
    # NOTE that using only_binary with this model will restrict much more
    # than using only_binary with the model without symmetry, unless the demand
    # of all pieces is naturally_only_binary the model will give much worse
    # results.
    @variable(model, pieces_sold[1:num_piece_types], Bin)
  else
    @variable(model,
      0 <= pieces_sold[pii = 1:num_piece_types] <= d[pii],
    Int)
  end

  if only_binary
    @variable(model, cuts_made[1:length(hvcuts)], Bin)
  else
    @variable(model, cuts_made[1:length(hvcuts)] >= 0, Int)
  end

  # The objective function maximizes the profit of the pieces sold.
  @objective(model, Max,
    sum(p[pii] * sum(pieces_sold[pii]) for pii = 1:num_piece_types)
  )

  # ub0: trivial area upper bound
  #@constraint(model, sum(a[j] * x[i, j] for i = 1:N, j = 1:T) <= L*W)

  # c1: There is just one of the original plate, and so it can be only used to
  # extract a single piece (this is rare, because there would need to be a
  # piece the exact same size as the original plate) xor make a single cut that
  # would make two new subplates available.
  @constraint(model,
    sum(pieces_sold[pli2piis[1]]) + sum(cuts_made[parent2cut[1]]) <= 1
  )

  # c2: for each subplate type that is not the original plate, such subplate
  # type may be cut at most the number of times it was a child of another cut.
  for pli in 2:num_plate_types
    @constraints model begin
      sum(cuts_made[parent2cut[pli]]) <= sum(cuts_made[child2cut[pli]])
    end
  end

  # c2.5: for each non-original plate type, the amount of times it can be cut
  # (in any way) is bounded by the amount of times such subplate that can be
  # cut from the original plate.
  # TODO: use the more intelligent bound that I have in checkvist to make
  # the constraint even tighter.
  # TODO: if some plate type leaves a considerable space at the borders
  # (i.e., replicating it for the bound leaves a loose bound) should we
  # add some other plates together in the border space just to strengthen
  # the bound (this would add more nonzeros).
  for ssplis in same_size_plis
    @assert !iszero(length(ssplis))
    @assert isone(length(unique!(map(i -> pli2lwsb[i][1], ssplis))))
    @assert isone(length(unique!(map(i -> pli2lwsb[i][2], ssplis))))
    @constraint(model,
      sum(cuts_made[ssplis]) <= pli2lwsb[ssplis[1]][4]
    )
  end

  # c3: finally, for each piece type, the amount of pieces sold of that type is
  # at most the number of plates with the piece exact size that were not cut to
  # make smaller plates.
  for pii in 1:num_piece_types
    @constraint(model, pieces_sold[pii] <= sum(cuts_made[vcat(child2cut[pii2plis[pii]]...)]) - sum(cuts_made[vcat(parent2cut[pii2plis[pii]]...)]))
  end

  model, hvcuts, pli2lwsb, pii2plis, pli2piis 
end # build_model_with_symmbreak

end # module
