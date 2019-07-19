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

  pli_lwb, hnnn, vnnn, np = partitions(P, d, l, w, L, W)
  nnn = vcat(hnnn, vnnn)

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

  num_plate_types = length(pli_lwb)

  # pli2pair: inverse index which given a plate index will return a list of all
  # picuts indexes (np variable) that cut some piece from some plate.
  pli2pair = [Vector{P}() for _ = 1:num_plate_types]
  # pii2pair: the same as pli2pair but for piece indexes.
  pii2pair = [Vector{P}() for _ = 1:num_piece_types]

  # The vectors below all have the same length as the number of plate types (to
  # allow indexing by plate type). The value of a position is a vector of
  # arbitrary length and irrelevant index, the values of this inner vector are
  # cut indexes. Such cut indexes are related to the plate that is the index of
  # the outer vector.
  if break_hvcut_symm
    # fchild_vcut[10], for example, has a vector of cut indexes in which the
    # plate type 10 appears as the first child of a vertical cut.
    # First CHILD of Horizontal CUT
    fchild_hcut = [Vector{P}() for _ = 1:num_plate_types]
    # First CHILD of Vertical CUT
    fchild_vcut = [Vector{P}() for _ = 1:num_plate_types]
    # Second CHILD of Any CUT
    schild_acut = [Vector{P}() for _ = 1:num_plate_types]
    # PARENT in an Horizontal CUT
    parent_hcut = [Vector{P}() for _ = 1:num_plate_types]
    # PARENT in a Vertical CUT
    parent_vcut = [Vector{P}() for _ = 1:num_plate_types]
    for i in eachindex(hnnn)
      n1, n2, n3 = hnnn[i]
      push!(parent_hcut[n1], i)
      push!(fchild_hcut[n2], i)
      push!(schild_acut[n3], i)
    end
    len_hnnn = length(hnnn)
    for i in eachindex(vnnn)
      n1, n2, n3 = vnnn[i]
      j = len_hnnn + i
      push!(parent_vcut[n1], j)
      push!(fchild_vcut[n2], j)
      push!(schild_acut[n3], j)
    end
  else
    # any CHILD plate to respective CUT indexes
    child2cut = [Vector{P}() for _ = 1:num_plate_types]
    # any PARENT plate to respective CUT indexes
    parent2cut = [Vector{P}() for _ = 1:num_plate_types]

    # Initialize all inverse indexes.
    for i in eachindex(nnn)
      n1, n2, n3 = nnn[i]
      push!(parent2cut[n1], i)
      push!(child2cut[n2], i)
      push!(child2cut[n3], i)
    end
  end
  for i in eachindex(np)
    pli, pii = np[i]
    push!(pli2pair[pli], i)
    push!(pii2pair[pii], i)
  end

  #@variable(model,
  #  0 <= hvcuts[i = 1:length(nnn)] <= pli_lwb[nnn[i][1]][3]
  #, Int)
  # If all pieces have demand one, a binary variable will suffice to make the
  # connection between a piece type and the plate it is extracted from.
  naturally_only_binary = all(di -> di <= 1, d)
  if naturally_only_binary
    if break_hvcut_symm
      @variable(model, picuts[1:length(np), 1:2], Bin)
    else
      @variable(model, picuts[1:length(np)], Bin)
    end
  elseif only_binary
    if break_hvcut_symm
      @variable(model, picuts[1:length(np), 1:2], Bin)
    else
      @variable(model, picuts[1:length(np)], Bin)
    end
  else
    # TODO: check if we can break the symmetry without doubling the number
    # of picut variables. This would need some change in the flow constraints.
    # Having two separate arbitrary variables for the same picut probably
    # creates the symmetries we are exactly trying to solve.
    if break_hvcut_symm
      @variable(model, picuts[1:length(np), 1:2] >= 0, Int)
    else
      @variable(model, picuts[1:length(np)] >= 0, Int)
    end
    #@variable(model,
    #  0 <= picuts[i = 1:length(np)] <= min(pli_lwb[np[i][1]][3], d[np[i][2]]),
    #Int)
  end

  if only_binary
    @variable(model, hvcuts[1:length(nnn)], Bin)
    if break_hvcut_symm
      @variable(model, hcuts[1:length(nnn)], Bin)
    end
  else
    @variable(model, hvcuts[1:length(nnn)] >= 0, Int)
    if break_hvcut_symm
      @variable(model, hcuts[1:length(nnn)] >= 0, Int)
    end
  end

  if break_hvcut_symm
    @objective(model, Max,
      sum(p[pii] * sum(picuts[pii2pair[pii], :]) for pii = 1:num_piece_types)
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
      sum(picuts[pli2pair[1], :]) +
      sum(hvcuts[vcat(parent_hcut[1], parent_vcut[1])]) == 1
    )
  else
    @constraint(model,
      sum(picuts[pli2pair[1]]) + sum(hvcuts[parent2cut[1]]) <= 1
    )
  end

  # c2: 
  for pli in 2:num_plate_types
    if break_hvcut_symm
      @constraints model begin
        sum(hcuts[schild_acut[pli]]) <= sum(hvcuts[schild_acut[pli]])
        sum(picuts[pli2pair[pli], 1]) + sum(hvcuts[parent_hcut[pli]]) <= (sum(hvcuts[fchild_vcut[pli]]) + sum(hvcuts[schild_acut[pli]]) - sum(hcuts[schild_acut[pli]]))
        sum(picuts[pli2pair[pli], 2]) + sum(hvcuts[parent_vcut[pli]]) <= (sum(hvcuts[fchild_hcut[pli]]) + sum(hcuts[schild_acut[pli]]))
      end
    else
      @constraint(model,
        sum(picuts[pli2pair[pli]]) + sum(hvcuts[parent2cut[pli]]) <=
        sum(hvcuts[child2cut[pli]])
      )
    end
  end

  # c2.5: The amount of each plate type generated by cuts is bounded by the
  # amount that can be cut from the original plate.
  for pli in 2:num_plate_types
    if break_hvcut_symm
      @constraint(model,
        sum(hvcuts[
          vcat(fchild_hcut[pli], fchild_vcut[pli], schild_acut[pli])
        ]) <= pli_lwb[pli][3]
      )
    else
      @constraint(model,
        sum(hvcuts[child2cut[pli]]) <= pli_lwb[pli][3]
      )
    end
  end

  # c3: 
  for pii in 1:num_piece_types
    if break_hvcut_symm
      @show pii
      @show pii2pair[pii]
      @constraint(model, sum(picuts[pii2pair[pii], :]) <= d[pii])
    else
      @constraint(model, sum(picuts[pii2pair[pii]]) <= d[pii])
    end
  end
  
  model, pli_lwb, nnn, np
end # build

end # module

