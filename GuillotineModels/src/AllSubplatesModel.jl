# Variables:
#
# `x[pli, pii, fci, lci]`: Integer. The number of subplates of type
# `pli` which have a piece of type `pii` cut from them and create the frontal
# and lateral childs `fci` and `lci`.
#
# Objective function:
#
# Maximize the profit of the pieces cut.
#   sum(p[pii] * x[_, pii, _, _])
#
#
# Constraints:
#
# There is just one of the original plate.
#   sum(x[1, _, _, _]) <= 1
# The number of subplates available depends on the number of plates that have
# it as children. Note: the below works because fci and pli are always
# different subplates.
#   sum(x[pli>1, _, _, _]) <= sum(x[_, _, fci, lci]) (fci or lci equals to pli)
# The number of pieces of some type is always less than or equal to the demand.
#   sum(x[_, pii, _, _]) <= d[pii]
#
#
# Unnecessary constraints:
#
# The number of times a pair pli-pii appear is at most the min between: 1)
# d[pii] and the number of subplates pli that fit in the original plate.
#   sum(x[pli>1, pii, _, _]) <= min(d[pii], max_fits[pli])
#
#
# The ideal pre-processing:
# 
# Given the constraints format to build the model with zero overhead we need:
#   An array of bounds 1:num_of_plate_types with arrays inside. The inner
#   arrays have a variable number of tuples (pii, fci, lci) in which
#   pii is unique (inside the inner vector).
#   The method will iterate this structure (inner first, outer second) and
#   for each element it will:
#     1) be added to the future lhs of constraint 1 or 2, the specific line
#        (and constraint family) will be given by pli.
#     2) be added to the future lhs of constraint 3, the specific line
#        will be given by pii;
#     3) be added to the future rhs of constraint 2, the specific lines will
#        be given by fci and lci;
#

module AllSubplatesModel

using JuMP
push!(LOAD_PATH, "./")
using GuillotinePlatesDP

# TODO: allow passing the partitions as parameters, decouple the generation
# of the inverse indexes and allow the the user to just generate the inverse
# indexes.
function build(
  model, d :: Vector{D}, p :: Vector{P}, l :: Vector{S}, w :: Vector{S},
  L :: S, W :: S
) where {D, S, P}
  @assert length(d) == length(l) && length(l) == length(w)
  num_piece_types = convert(D, length(d))

  num_plate_types, parts = partitions(D, P, l, w, L, W)
  num_parts = length(parts)

  # A Vector of the same length as the number of plate types, in which the values
  # are Vectors of partition indexes in which the plate cut is the outer Vector
  # index.
  c1and2lhss = [Vector{P}() for _ = 1:num_plate_types]
  # A Vector of the same length as the number of plate types, in which the values
  # are Vectors of partition indexes in which one of the child plates is the
  # outer Vector index.
  c2rhss = [Vector{P}() for _ = 1:num_plate_types]
  # A Vector of length equal to the number of piece types, in which the values
  # are Vectors of the indexes of the partitions that cut such piece type.
  c3lhss = [Vector{P}() for _ = 1:num_piece_types]

  # Fills the three structures below. The three of them are just inverted indexes.
  for pai in 1:num_parts
    pli, pii, fci, lci = parts[pai]
    push!(c1and2lhss[pli], pai)
    push!(c3lhss[pii], pai)
    !iszero(fci) && push!(c2rhss[fci], pai)
    !iszero(lci) && push!(c2rhss[lci], pai)
  end
  # @assert explanation: No partition can have a child plate with the
  # dimensions of the original plate.
  @assert iszero(length(c2rhss[1]))

  @variable(model, x[1:num_parts], Int)

  @objective(model, Max,
    sum(p[pii] * sum(x[c3lhss[pii]]) for pii = 1:num_piece_types)
  )

  # ub0: trivial area upper bound
  #@constraint(model, sum(a[j] * x[i, j] for i = 1:N, j = 1:T) <= L*W)

  # c1: 
  @constraint(model, sum(x[c1and2lhss[1]]) <= 1)

  # c2: 
  for pli in 2:num_plate_types
    @constraint(model, sum(x[c1and2lhss[pli]]) <= sum(x[c2rhss[pli]]))
  end

  # c3: 
  for pii in 1:num_piece_types
    @constraint(model, sum(x[c3lhss[pii]]) <= d[pii])
  end
  
  model, x
end # build

end # module

