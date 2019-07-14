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
  L :: S, W :: S
) where {D, S, P}
  @assert length(d) == length(l) && length(l) == length(w)
  num_piece_types = convert(D, length(d))

  ilwb, hnnn, vnnn, np = partitions(P, d, l, w, L, W)
  num_plate_types = length(ilwb)
  ilwb_ = Vector{Tuple{S, S, D}}(undef, num_plate_types)
  for (pli, pll, plw, plb) in ilwb
    ilwb_[pli] = (pll, plw, plb)
  end

  # A Vector of the same length as the number of plate types, in which the values
  # are Vectors of partition indexes in which the plate cut is the outer Vector
  # index.
  #=
  pli2hpair = [Vector{P}() for _ = 1:num_plate_types]
  pli2vpair = [Vector{P}() for _ = 1:num_plate_types]
  =#
  pli2pair = [Vector{P}() for _ = 1:num_plate_types]
  # 
  # The vectors below all have the same length as the number of plate types (
  # to allow indexing by plate type). The value of a position is a vector of
  # arbitrary length and irrelevant index and its values are cut indexes.
  # The cut indexes of the inner vector are related to the plate that is the
  # index of the outer vector with the inner vector inside.
  # fchild_vcut[10], for example, has a vector of cut indexes in which the
  # plate type 10 appears as the first child of the cut
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
  # A Vector of length equal to the number of piece types, in which the values
  # are Vectors of the indexes of the partitions that cut such piece type.
  pii2pair = [Vector{P}() for _ = 1:num_piece_types]

  # Fills the inverted indexes.
  for i in eachindex(np)
    pli, pii = np[i]
    #=
    ilwb_[pli][1] < l[pii] && @show pli, pii, ilwb_[pli][1], l[pii]
    ilwb_[pli][2] < w[pii] && @show pli, pii, ilwb_[pli][2], w[pii]
    if ilwb_[pli][1] - l[pii] >= ilwb_[pli][2] - w[pii] 
      push!(pli2vpair[pli], i)
    else
      push!(pli2hpair[pli], i)
    end
    =#
    push!(pli2pair[pli], i)
    push!(pii2pair[pii], i)
  end

  for i in eachindex(hnnn)
    n1, n2, n3 = hnnn[i]
    push!(parent_hcut[n1], i)
    push!(fchild_vcut[n2], i)
    push!(schild_acut[n3], i)
  end
  len_hnnn = length(hnnn)
  for i in eachindex(vnnn)
    n1, n2, n3 = vnnn[i]
    j = len_hnnn + i
    push!(parent_vcut[n1], j)
    push!(fchild_hcut[n2], j)
    push!(schild_acut[n3], j)
  end
  nnn = vcat(hnnn, vnnn)

  # If all pieces have demand one, a binary variable will suffice to make the
  # connection between a piece type and the plate it is extracted from.
  if all(di -> di <= 1, d)
    @variable(model, picuts[1:length(np), 1:2], Bin)
  else
    @variable(model, picuts[1:length(np), 1:2] >= 0, Int)
  end
  @variable(model, hvcuts[1:length(nnn)] >= 0, Int)
  @variable(model, hcuts[1:length(nnn)] >= 0, Int)

  @objective(model, Max,
    sum(p[pii] * sum(picuts[pii2pair[pii], :]) for pii = 1:num_piece_types)
  )

  # ub0: trivial area upper bound
  #@constraint(model, sum(a[j] * x[i, j] for i = 1:N, j = 1:T) <= L*W)

  # c1: 
  @constraint(model, sum(picuts[pli2pair[1], :]) + sum(hvcuts[vcat(parent_hcut[1], parent_vcut[1])]) == 1)

  # c2: 
  for pli in 2:num_plate_types
    @constraints model begin
      sum(hcuts[schild_acut[pli]]) <= sum(hvcuts[schild_acut[pli]])
      sum(picuts[pli2pair[pli], 1]) + sum(hvcuts[parent_hcut[pli]]) <= (sum(hvcuts[fchild_vcut[pli]]) + sum(hvcuts[schild_acut[pli]]) - sum(hcuts[schild_acut[pli]]))
      sum(picuts[pli2pair[pli], 2]) + sum(hvcuts[parent_vcut[pli]]) <= (sum(hvcuts[fchild_hcut[pli]]) + sum(hcuts[schild_acut[pli]]))
    end
  end

  # c2.5: The amount of each plate type generated by cuts is bounded by the
  # amount that can be cut from the original plate.
  #for pli in 2:num_plate_types
  #  @constraint(model, sum(hvcuts[vcat(fchild_vcut[pli], fchild_hcut[pli], schild_acut[pli])]) <= ilwb_[pli][3])
  #end

  # c3: 
  for pii in 1:num_piece_types
    @constraint(model, sum(picuts[pii2pair[pii], :]) <= d[pii])
  end
  
  model, ilwb, nnn, np
end # build

end # module

