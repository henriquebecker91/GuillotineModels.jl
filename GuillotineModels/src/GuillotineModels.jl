module GuillotineModels

using JuMP

# d: piece demand (d)
# l: piece dimension (l or w)
# L: plate dimension (L or W)
function becker2019_discretize(
  d :: Vector{D}, l :: Vector{S}, L :: S
) where {D, S}
  # If two pieces have the same dimension they should be merged into a single
  # piece (summing their demands) for performance reasons.
  @assert l == unique(l)
  # Each piece in 1:N has a demand and a dimension.
  @assert length(d) == length(l)
  N = length(d)
  # marks: for each unit of the plate dimension, if there can be a cut there 
  # (considering the pieces available) or not.
  marks = fill(false, L)
  # Mark the cuts of the first piece.
  marks[l[1]] = true
  y = l[1] # y: used to iterate capacity, inherited from knapsack papers
  for _ = 2:d[1]
    marks[y += l[1]] = true
  end
  # Mark the cuts of all other pieces.
  for pii = 2:N # PIece Index
    li = l[pii] # length of i
    di = d[pii] # demand of i
    for y = (L - li):-1:1
      if marks[y]
        yrli = y + li # yrli: y + repeated lengths of i
        marks[yrli] = true
        for r = 2:di
          yrli += li
          yrli > L && break
          marks[yrli] = true
        end
      end
    end
    # Mark cuts considering the pieces began to (could be done using a dummy cut at
    # index zero with value true, but the array has no index zero). 
    y = li
    marks[y] = true
    for _ = 2:di
      marks[y += li] = true
    end
  end

  cuts = Vector{S}()
  sizehint!(cuts, L)
  for (position, is_marked) in enumerate(marks)
    is_marked && push!(cuts, position)
  end
  cuts
end

function herz_discretize(d, D)
  P = 
  c = zeros(size(D+1))

  for i = 1:length(d) 
    for j = d[i]:D 
      if c[j+1] < c[j+1-d[i]] + d[i] 
        c[j+1] = c[j+1-d[i]] + d[i]
      end
    end
  end

  for j = 1:D 
    if c[j+1] == j 
      push!(P, j)
    end
  end

  c
end

function mock_discretize(d :: Vector{T}, D :: T) where {T}
  return collect(one(T):D), D
end

function build(model, l, w, d, L, W; p = nothing) 
  @assert length(l) == length(w) && length(w) == length(d)
  T = length(l)
  N = 2*T + 1
  a = l .* w # area of the piece types
  dl, DL = discretize(l, L)
  dw, DW = discretize(w, W)

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
