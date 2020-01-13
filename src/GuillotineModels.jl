module GuillotineModels

using JuMP

function build(
  model, ::Type{P}, d :: Vector{D}, l :: Vector{S}, w :: Vector{S},
  L :: S, W :: S
) where {D, S, P}
  @assert length(l) == length(w) && length(w) == length(d)
  a = l .* w # area of the piece types
  
  partitions(P, d, l, w, L, W)
  @variables model begin
    x[1:N, 1:T], Bin  # true if piece type is placed at node, false otherwise
    v[1:N, 1:DL], Bin # true if node is a vertical cut, false otherwise
    h[1:N, 1:DW], Bin # true if node is a horizontal cut, false otherwise
    l_[1:N] >= 0 # length available to node N
    w_[1:N] >= 0 # width available to node N
  end

  @objective(model, Max, sum(p[j] * x[i, j] for i = 1:N, j = 1:T))

  # ub0: trivial area upper bound
  @constraint(model, sum(a[j] * x[i, j] for i = 1:N, j = 1:T) <= L*W)

  # c1: Each node is either: 1) a placed piece (leaf); 2) a horizontal cut
  # (branch); 3) a vertical cut (branch).
  for i = 1:N
    @constraint(model,
      sum(x[i, j] for j = 1:T) + sum(v[i, j] for j = 1:DL) +
      sum(h[i, j] for j = 1:DW) <= 1
    )
  end
  # c2: Each non-root node can only be a branch or leaf if the node parent
  # is a branch.
  for i = 2:N
    @constraint(model,
      sum(x[i, j] for j = 1:T) + sum(v[i, j] for j = 1:DL) +
      sum(h[i, j] for j = 1:DW) <= sum(v[div(i, 2), j] for j = 1:DL) +
      sum(h[div(i, 2), j] for j = 1:DW)
    )
  end

  # c3 and c4: The root node is limited by the size of the original plate.
  fix(l_[1], L; force = true)
  fix(w_[1], W; force = true)

  # c5: There is a limit on the number of pieces available for each piece type.
  for j = 1:T
    @constraint(model,
      sum(x[i, j] for i = 1:N) <= d[j]
    )
  end

  # c6 and c7: The cuts reduce the space available for their children.
  # c9 and c10: the length (width) of a plate is always at max the length
  # NOTE: it is not possible to tighten the big-M using dl[1] and dw[1] because
  # they can be children that yet keep one dimension at the same size as root.
  for i = 2:N
    i2 = div(i, 2)
    if iseven(i) # the first child needs a big-M to restrict its dimensions
      @constraints model begin
        l_[i] <= sum(dl[j] * v[i2, j] for j = 1:DL) + L * sum(h[i2, j] for j = 1:DW)
        w_[i] <= sum(dw[j] * h[i2, j] for j = 1:DW) + W * sum(v[i2, j] for j = 1:DL)
        # these are c9 and c10, they are implicit for the isodd(i) case
        l_[i] <= l_[i2]
        w_[i] <= w_[i2]
      end
    else
      # the second child does not need a big-M to restrict its dimensions,
      # nor does it need c9 and c10
      @constraints model begin
        l_[i] <= l_[i2] - sum(dl[j] * v[i2, j] for j = 1:DL)
        w_[i] <= w_[i2] - sum(dw[j] * h[i2, j] for j = 1:DW)
      end
    end
  end
  # c11 and c12: the placed pieces respect length (width) of the plate where
  # they are placed.
  for i = 1:N
    @constraints model begin
      sum(l[j] * x[i, j] for j = 1:T) <= l_[i]
      sum(w[j] * x[i, j] for j = 1:T) <= w_[i]
    end
  end
  
  model
end # build

end # module
