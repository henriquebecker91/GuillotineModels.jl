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

  ilwb, nnn, np = partitions(P, d, l, w, L, W)
  num_plate_types = length(ilwb)

  # A Vector of the same length as the number of plate types, in which the values
  # are Vectors of partition indexes in which the plate cut is the outer Vector
  # index.
  pli2pair = [Vector{P}() for _ = 1:num_plate_types]
  # A Vector of the same length as the number of plate types, in which the values
  # are Vectors of partition indexes in which one of the child plates is the
  # outer Vector index.
  child2cut = [Vector{P}() for _ = 1:num_plate_types]
  parent2cut = [Vector{P}() for _ = 1:num_plate_types]
  # A Vector of length equal to the number of piece types, in which the values
  # are Vectors of the indexes of the partitions that cut such piece type.
  pii2pair = [Vector{P}() for _ = 1:num_piece_types]

  # Fills the three structures below.
  # The three of them are just inverted indexes.
  for i in eachindex(nnn)
    n1, n2, n3 = nnn[i]
    push!(parent2cut[n1], i)
    push!(child2cut[n2], i)
    push!(child2cut[n3], i)
  end
  for i in eachindex(np)
    pli, pii = np[i]
    push!(pli2pair[pli], i)
    push!(pii2pair[pii], i)
  end

  # @assert explanation: No cut can generate a child plate with the
  # dimensions of the original plate.
  @assert iszero(length(child2cut[1]))

  @variables model begin
    picuts[1:length(np)] >= 0, Int
    hvcuts[1:length(nnn)] >= 0, Int
  end

  @objective(model, Max,
    sum(p[pii] * sum(picuts[pii2pair[pii]]) for pii = 1:num_piece_types)
  )

  # ub0: trivial area upper bound
  #@constraint(model, sum(a[j] * x[i, j] for i = 1:N, j = 1:T) <= L*W)

  # c1: 
  @constraint(model, sum(picuts[pli2pair[1]]) + sum(hvcuts[parent2cut[1]]) == 1)

  # c2: 
  for pli in 2:num_plate_types
    @constraint(model, sum(picuts[pli2pair[pli]]) + sum(hvcuts[parent2cut[pli]]) <= sum(hvcuts[child2cut[pli]]))
  end

  # c3: 
  for pii in 1:num_piece_types
    @constraint(model, sum(picuts[pii2pair[pii]]) <= d[pii])
  end
  
  model, ilwb, nnn, np
end # build

end # module

