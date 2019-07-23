# Variables:
#
# `hvcuts[n1, n2, n3]`: Integer. The number of subplates of type
# `n1` that are cut into subplates `n2` and `n3` (horizontal and
# vertical cuts are together for now).
# `picuts[n, p]`: Integer. The number of subplates of type `n` that
# were used to generate piece `p`.
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
#   sum(picuts[1, _]) + sum(hvcuts[1, _,  _]) = 1
# The number of subplates available depends on the number of plates that have
# it as children.
#   sum(picuts[n1>1, _]) + sum(hvcuts[n1>1, _, _]) <= sum(hvcuts[_, n2, n3])
#     where n2 == n1 or n3 == n1, double hvcuts[_, n2, n3] if n2 == n3 
# The number of pieces of some type is always less than or equal to the demand.
#   sum(picuts[_, pii]) <= d[pii]
#
# Unnecessary constraints:
#
# The number of times a pair pli-pii appear is at most the min between: 1)
# d[pii] and the number of subplates pli that fit in the original plate.
#   sum(picuts[n, pii]) <= min(d[pii], max_fits[n])

module AllSubplatesModel

using JuMP
push!(LOAD_PATH, "./")
using GuillotinePlatesDP

function build(
  model, d :: Vector{D}, p :: Vector{P}, l :: Vector{S}, w :: Vector{S},
  L :: S, W :: S; only_binary = false, break_hvcut_symm = false
) where {D, S, P}
  @assert length(d) == length(l) && length(l) == length(w)
  num_piece_types = convert(D, length(d))

  sllw = SortedLinkedLW(D, l, w)
  if break_hvcut_symm
    pli2lwsb, hcuts, vcuts, pii2plis, pli2piis = partitions_no_symm(
      P, d, sllw, L, W
    )
    num_plate_types = length(pli2lwsb)
  else
    pli_lwb, hcuts, vcuts, np = partitions(P, d, sllw, L, W)
    num_plate_types = length(pli_lwb)
  end
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

  if !break_hvcut_symm
    # pli2pair: inverse index which given a plate index will return a list of
    # all picuts indexes (np variable) that cut some piece from some plate.
    pli2pair = [Vector{P}() for _ = 1:num_plate_types]
    # pii2pair: the same as pli2pair but for piece indexes.
    pii2pair = [Vector{P}() for _ = 1:num_piece_types]
  end

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

  if !break_hvcut_symm
    for i in eachindex(np)
      pli, pii = np[i]
      push!(pli2pair[pli], i)
      push!(pii2pair[pii], i)
    end
  end

  #@variable(model,
  #  0 <= hvcuts[i = 1:length(nnn)] <= pli_lwb[nnn[i][1]][3]
  #, Int)
  # If all pieces have demand one, a binary variable will suffice to make the
  # connection between a piece type and the plate it is extracted from.
  naturally_only_binary = all(di -> di <= 1, d)
  if naturally_only_binary
    if break_hvcut_symm
      @variable(model, pieces_sold[1:num_piece_types], Bin)
    else
      @variable(model, picuts[1:length(np)], Bin)
    end
  elseif only_binary # is the same as naturally_only_binary, but only for now
    if break_hvcut_symm
      @variable(model, pieces_sold[1:num_piece_types], Bin)
    else
      @variable(model, picuts[1:length(np)], Bin)
    end
  else
    if break_hvcut_symm
      @variable(model,
        0 <= pieces_sold[pii = 1:num_piece_types] <= d[pii],
      Int)
    else
      @variable(model, picuts[1:length(np)] >= 0, Int)
      #@variable(model,
      #  0 <= picuts[i = 1:length(np)] <=
      #    min(pli_lwb[np[i][1]][3], d[np[i][2]]),
      #Int)
    end
  end

  if only_binary
    @variable(model, cuts_made[1:length(hvcuts)], Bin)
  else
    @variable(model, cuts_made[1:length(hvcuts)] >= 0, Int)
  end

  if break_hvcut_symm
    @objective(model, Max,
      sum(p[pii] * sum(pieces_sold[pii]) for pii = 1:num_piece_types)
    )
  else
    @objective(model, Max,
      sum(p[pii] * sum(picuts[pii2pair[pii]]) for pii = 1:num_piece_types)
    )
  end

  # ub0: trivial area upper bound
  #@constraint(model, sum(a[j] * x[i, j] for i = 1:N, j = 1:T) <= L*W)

  # c1: 
  if break_hvcut_symm
    @constraint(model,
      sum(pieces_sold[pli2piis[1]]) + sum(cuts_made[parent2cut[1]]) <= 1
    )
  else
    @constraint(model,
      sum(picuts[pli2pair[1]]) + sum(cuts_made[parent2cut[1]]) <= 1
    )
  end

  # c2: 
  for pli in 2:num_plate_types
    if break_hvcut_symm
      @constraints model begin
        sum(cuts_made[parent2cut[pli]]) <= sum(cuts_made[child2cut[pli]])
      end
    else
      @constraint(model,
        sum(picuts[pli2pair[pli]]) + sum(cuts_made[parent2cut[pli]]) <=
        sum(cuts_made[child2cut[pli]])
      )
    end
  end

  # c2.5: The amount of each plate type generated by cuts is bounded by the
  # amount that can be cut from the original plate.
  # TODO: do this in a more intelligent fashion, group the possibly three
  # plate codes that have the same exact size (reduce the number of
  # constraints by a factor of three and make them tighter).
  # TODO: use the more intelligent bound that I have in checkvist to make
  # the constraint even tighter.
  # TODO: if some plate type leaves a considerable space at the borders
  # (i.e., replicating it for the bound leaves a loose bound) should we
  # add some other plates together in the border space just to strengthen
  # the bound (this would add more nonzeros).
  for pli in 2:num_plate_types
    if break_hvcut_symm
      @constraint(model,
        sum(cuts_made[parent2cut[pli]]) <= pli2lwsb[pli][4]
      )
    else
      @constraint(model,
        sum(cuts_made[parent2cut[pli]]) <= pli_lwb[pli][3]
      )
    end
  end

  # c3: 
  for pii in 1:num_piece_types
    if break_hvcut_symm
      @constraint(model, pieces_sold[pii] <= sum(cuts_made[vcat(child2cut[pii2plis[pii]]...)]) - sum(cuts_made[vcat(parent2cut[pii2plis[pii]]...)]))
    else
      @constraint(model, sum(picuts[pii2pair[pii]]) <= d[pii])
    end
  end
  
  #model, hvcuts, pli_lwb, np
  model, hvcuts, pli2lwsb
end # build

end # module

