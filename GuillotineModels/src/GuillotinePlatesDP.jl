module GuillotinePlatesDP

push!(LOAD_PATH, "./")
using GC2DInstanceReader

export all_plates, nibnn, write_nibnn, partitions

function fits(l :: S, w :: S, L :: S, W :: S) :: Bool where {S}
  l <= L && w <= W
end

# d: piece demand (d)
# l: piece dimension (l or w)
# L: plate dimension (L or W)
function becker2019_discretize(
  d :: Vector{D}, l :: Vector{S}, L :: S
) where {D, S}
  # If two pieces have the same dimension they should be merged into a single
  # piece (summing their demands) for performance reasons.
  #@assert l == unique(l)
  # Each piece in 1:N has a demand and a dimension.
  @assert length(d) == length(l)
  N = length(d)
  # Remove items that cannot fit into L anyway.
  d_, l_ = similar(d), similar(l)
  j = 0
  for i in 1:N
    if l[i] <= L
      j += 1
      d_[j] = d[i]
      l_[j] = l[i]
    end
  end
  resize!(d_, j)
  resize!(l_, j)
  d, l = d_, l_
  N = j
  iszero(N) && return S[]
  # marks: for each unit of the plate dimension, if there can be a cut there 
  # (considering the pieces available) or not.
  marks = fill(false, L)
  # Mark the cuts of the first piece.
  marks[l[1]] = true
  y = l[1] # y: used to iterate capacity, inherited from knapsack papers
  for _ = 2:min(d[1], L ÷ l[1])
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
    for _ = 2:min(di, L ÷ li)
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

function becker2019_discretize(
  d :: Vector{D}, l :: Vector{S}, w :: Vector{S}, L :: S, W :: S
) where {D, S}
  # If two pieces have the same dimension they should be merged into a single
  # piece (summing their demands) for performance reasons.
  #@assert l == unique(l)
  # Each piece in 1:N has a demand and a dimension.
  @assert length(d) == length(l)
  N = length(d)
  # Remove items that cannot fit into L anyway.
  d_, l_, w_ = similar.((d, l, w))
  j = 0
  for i in 1:N
    if l[i] <= L && w[i] <= W
      j += 1
      d_[j] = d[i]
      l_[j] = l[i]
      w_[j] = w[i]
    end
  end
  resize!(d_, j)
  resize!(l_, j)
  resize!(w_, j)
  d, l, w = d_, l_, w_
  N = j
  iszero(N) && return S[]
  # marks: for each unit of the plate dimension, if there can be a cut there 
  # (considering the pieces available) or not.
  marks = fill(false, L)
  # Mark the cuts of the first piece.
  marks[l[1]] = true
  y = l[1] # y: used to iterate capacity, inherited from knapsack papers
  for _ = 2:min(d[1], L ÷ l[1])
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
    for _ = 2:min(di, L ÷ li)
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

function is_leaf_plate(
  l :: S, w :: S, L :: S, W :: S, minl :: S, minw :: S
) :: Bool where {S}
  L >= l && W >= w && l < (L + minl) && w < (W + minw) 
end

function partitions(
  ::Type{P}, d :: Vector{D}, l :: Vector{S}, w :: Vector{S}, L :: S, W :: S
) where {D, S, P}
  @assert length(d) == length(l)
  @assert length(d) == length(w)
  max_piece_type = convert(D, length(l))
  max_num_plates = convert(P, L) * convert(P, W)
  #dl = becker2019_discretize(d, l, L ÷ 2)
  #dw = becker2019_discretize(d, w, W ÷ 2)
  min_pil = minimum(l)
  min_piw = minimum(w)
  # If (n1, n2, n3) is in nnn, then n1 may be partitioned in n2 and n3.
  # Not every possible partition is present, just the ones which
  # each child plate can fit at least one piece.
  hnnn = Vector{NTuple{3, P}}()
  vnnn = Vector{NTuple{3, P}}()
  # If the piece fits the plate size and the plate does not allow any other
  # piece by its side (i.e., there is no second piece that can be placed in the
  # same plate vertically or horizontally) then (plate, piece) is in np.
  np = Vector{Tuple{P, D}}()
  # The list of plates attributes: plate index, plate length, plate width,
  # and plate bound.
  ilwb = Vector{Tuple{P, S, S, P}}()
  # plis: matrix of the plate dimensions in which zero means "never seen that
    # plate before" and nonzero means "this nonzero number is the plate index".
  plis = zeros(P, L, W)
  plis[L, W] = one(P)
  # next: plates already indexed but not yet processed, starts with (L, W, 1).
  # Storing the index as the third value is not necessary (as it could be
  # queried from plis) but this is probably more efficient this way.
  next = Vector{Tuple{S, S, P}}()
  #sizehint!(next, max_num_plates)
  push!(next, (L, W, one(P)))
  # n: The amount of plates (the index of the highest plate type).
  n = one(P) # there is already the original plate
  while !isempty(next)
    pll, plw, pli = pop!(next) # PLate Length, Width, and Index
    plb = (L ÷ pll) * (W ÷ plw) # PLate Bound
    push!(ilwb, (pli, pll, plw, plb))
    for pii in 1:max_piece_type # pii: PIece Index
      pil, piw = l[pii], w[pii] # PIece Length and Width (homophones, I know)
      if is_leaf_plate(pil, piw, pll, plw, min_pil, min_piw)
        push!(np, (pli, pii))
      end
    end
    #for y in dl
    for y in becker2019_discretize(d, l, w, pll ÷ 2, plw)
      #y > pll ÷ 2 && break
      @assert plw >= min_piw
      @assert y >= min_pil
      @assert pll - y >= min_pil
      if iszero(plis[y, plw])
        push!(next, (y, plw, n += 1))
        plis[y, plw] = n
      end
      if iszero(plis[pll - y, plw])
        push!(next, (pll - y, plw, n += 1))
        plis[pll - y, plw] = n
      end
      push!(hnnn, (pli, plis[y, plw], plis[pll - y, plw]))
    end
    #for x in dw
    for x in becker2019_discretize(d, w, l, plw ÷ 2, pll)
      #x > plw ÷ 2 && break
      @assert pll >= min_pil
      @assert x >= min_piw
      @assert plw - x >= min_piw
      if iszero(plis[pll, x])
        push!(next, (pll, x, n += 1))
        plis[pll, x] = n
      end
      if iszero(plis[pll, plw - x])
        push!(next, (pll, plw - x, n += 1))
        plis[pll, plw - x] = n
      end
      push!(vnnn, (pli, plis[pll, x], plis[pll, plw - x]))
    end
  end

  return ilwb, hnnn, vnnn, np
end

end # module

