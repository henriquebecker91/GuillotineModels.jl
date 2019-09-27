module AllSubplatesModel

using JuMP
include("GuillotinePlatesDP.jl")
using .GuillotinePlatesDP

function search_cut(
  pp :: P, # parent plate
  fcl :: S, # first child length
  fcw :: S, # first child width
  nnn :: Vector{NTuple{3, P}},
  pli_lwb :: Vector{Tuple{S, S, P}},
) :: P where {D, S, P}
  for i in 1:length(nnn)
    (pp_, fc, _) = nnn[i]
    pp_ != pp && continue
    pli_lwb[fc][1] == fcl && pli_lwb[fc][2] == fcw && return 
  end
  @assert false # this should not be reachable
end

# TODO: the warm-start for the flag faithful2furini enabled and disabled will
# need to be different? the rules for which plates exist and which do not are
# different from one to another.
# NOTE: this method only work for simple patterns in which:
# (i) the cuts are two-staged (i.e., the pattern is justa a vector of vector);
# (ii) the first stage cuts vertically (width strips);
# (iii) the first piece of each strip gives the width of the strip;
# If you need to warm-start with a more complex pattern, create another
# method with the same name, and another type for parameter `pat`.
function warm_start(
  model, l, w, L, W,
  pat :: AbstractVector{AbstractVector{D}},
  pli_lwb :: Vector{Tuple{S, S, P}},
  nnn :: Vector{NTuple{3, P}},
  np :: Vector{Tuple{P, D}},
  plis :: Array{P, 2},
  faithful2furini = false
  #round2disc wait to see if this is needed
  # which other model building options will need to be passed to this?
) :: where {D, S, P}
  # the initial residual plate is L, W
  # visit the outer vector in reverse
  # if the current head of stripe is smaller than half residual plate
  # then search for a vertical cut on the residual plate, with the right width
  #   for the first child and enable it, change the residual plate to
  #   be the second child
  # else assert this is the last stripe, just use the remaining plate (second
  #   child of the last cut, or the whole root if there is just one stripe)
  # after finishing the strip processing, for each plate that is a stripe:
  #   do the same as the first stage, but for the subplate and the opposite cut
  #   orientation (including iteration in reverse order); 
  # finally, for every subplate that will become a piece, connect it to a piece
  #   for faithful2furini2016 we need to trim the plate and have it with exact
  #   plate size; for !faithful2furini2016 we just lik directly to np
  rl, rw = L, W
  rpli = plis[rl, rw]
  final = false
  for stripe in reverse(pat)
    @assert !final
    @assert !isempty(stripe)
    ws = w[first(stripe)]
    if ws > div(rw, 2)
      final = true
    else
      
    end
  end
end

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
  L :: S, W :: S; only_binary = false, use_c25 = false,
  ignore_2th_dim = false, ignore_d = false, round2disc = true,
  faithful2furini2016 = false,
  no_redundant_cut = false, no_cut_position = false,
  no_furini_symmbreak = false,
  relax2lp = false,
  lb :: P = zero(P), ub :: P = zero(P)
) where {D, S, P}
  @assert length(d) == length(l) && length(l) == length(w)
  num_piece_types = convert(D, length(d))

  sllw = SortedLinkedLW(D, l, w)
  pli_lwb, hcuts, vcuts, np = gen_cuts(P, d, sllw, L, W;
    ignore_2th_dim = ignore_2th_dim,
    ignore_d = ignore_d,
    round2disc = round2disc,
    no_cut_position = no_cut_position,
    no_redundant_cut = no_redundant_cut,
    no_furini_symmbreak = no_furini_symmbreak,
    faithful2furini2016 = faithful2furini2016
  )
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
    #=if iszero(schild)
      @show parent
      @show pli_lwb[parent]
      @show fchild
      @show pli_lwb[fchild]
      @show schild
    end=#
    #@assert faithful2furini2016 || !iszero(schild)
    !iszero(schild) && push!(child2cut[schild], i)
  end

  for i in eachindex(np)
    pli, pii = np[i]
    push!(pli2pair[pli], i)
    push!(pii2pair[pii], i)
  end

  # If all pieces have demand one, a binary variable will suffice to make the
  # connection between a piece type and the plate it is extracted from.
  naturally_only_binary = all(di -> di <= 1, d)
  if !relax2lp && (naturally_only_binary || only_binary)
    # only_binary is equal to naturally_only_binary for now, but the idea is
    # that only_binary will expand the number of binary variables to account
    # for picuts that can repeat (plates that can appear more than one time
    # and that may have the same piece extracted from them)
    @variable(model, picuts[1:length(np)], Bin)
  else
    @variable(model, picuts[1:length(np)] >= 0, integer = !relax2lp)
    #@variable(model,
    #  0 <= picuts[i = 1:length(np)] <=
    #    min(pli_lwb[np[i][1]][3], d[np[i][2]]),
    #Int)
  end

  if only_binary
    @assert !relax2lp
    @variable(model, cuts_made[1:length(hvcuts)] >= 0, binary = !relax2lp)
  else
    @variable(model, cuts_made[1:length(hvcuts)] >= 0, integer = !relax2lp)
  end

  # The objective function is to maximize the profit made by extracting
  # pieces from subplates.
  @objective(model, Max,
    sum(p[pii] * sum(picuts[pii2pair[pii]]) for pii = 1:num_piece_types)
  )

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

  if use_c25
    # c2.5: The amount of each subplate type generated by cuts (and used either
    # as a piece or as a intermediary plate) is bounded by the amount that can be
    # cut from the original plate.
    for pli in 2:num_plate_types
      @constraint(model,
        sum(picuts[pli2pair[pli]]) + sum(cuts_made[parent2cut[pli]]) <=
        pli_lwb[pli][3]
      )
    end
  end

  # c3: the amount of each piece type extracted from different plate types
  # cannot surpass the demand for that piece type.
  for pii in 1:num_piece_types
    @constraint(model, sum(picuts[pii2pair[pii]]) <= d[pii])
  end

  if !iszero(lb)
    @constraint(model,
      sum(p[pii]*sum(picuts[pii2pair[pii]]) for pii = 1:num_piece_types) >= (lb + 1)
    )
  end

  if !iszero(ub) && ub < sum(d .* p)
    @constraint(model,
      sum(p[pii]*sum(picuts[pii2pair[pii]]) for pii = 1:num_piece_types) <= ub
    )
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
  L :: S, W :: S; only_binary = false, use_c25 = false,
  ignore_2th_dim = false, ignore_d = false, round2disc = true
) where {D, S, P}
  @assert length(d) == length(l) && length(l) == length(w)
  num_piece_types = convert(D, length(d))

  sllw = SortedLinkedLW(D, l, w)
  pli2lwsb, hcuts, vcuts, pii2plis, pli2piis, same_size_plis =
    gen_cuts_sb(P, d, sllw, L, W; ignore_2th_dim = ignore_2th_dim,
    ignore_d = ignore_d,
    round2disc = round2disc
  )
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

  if use_c25
    @error "not sure if correctly implemented, check before"
    # TODO: check if the below is corret. What seems to be wrong is that ssplis
    # is a vector of plate indexes, while cuts_made should be indexed by
    # vectors of cut indexes.

    # c2.5: The amount of each subplate type generated by cuts (and used either
    # as a piece or as a intermediary plate) is bounded by the amount that can be
    # cut from the original plate.
    for ssplis in same_size_plis
      @assert !iszero(length(ssplis))
      @assert isone(length(unique!(map(i -> pli2lwsb[i][1], ssplis))))
      @assert isone(length(unique!(map(i -> pli2lwsb[i][2], ssplis))))
      @constraint(model,
        sum(cuts_made[ssplis]) <= pli2lwsb[ssplis[1]][4]
      )
    end
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

