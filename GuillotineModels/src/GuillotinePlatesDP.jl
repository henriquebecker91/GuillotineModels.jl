module GuillotinePlatesDP

push!(LOAD_PATH, "./")
using GC2DInstanceReader

export all_plates, nibnn, write_nibnn, partitions, partitions_no_symm, SortedLinkedLW

struct SortedLinkedLW{D, S}
  l :: Vector{S}
  w :: Vector{S}
  sl :: Vector{S}
  sw :: Vector{S}
  sli2pii :: Vector{D} # sorted length index to piece index
  swi2pii :: Vector{D} # sorted width index to piece index
  pii2sli :: Vector{D} # piece index to sorted length index
  pii2swi :: Vector{D} # piece index to sorted width index
end

function SortedLinkedLW(::Type{D}, l :: Vector{S}, w :: Vector{S}) where {D, S}
  @assert length(l) == length(w)
  n = length(l)
  sl = sort(l)
  sw = sort(w)
  sli2pii = sort!(collect(1:n), by = pii -> l[pii])
  swi2pii = sort!(collect(1:n), by = pii -> w[pii])
  pii2sli = Vector{D}(undef, n)
  pii2swi = Vector{D}(undef, n)
  for si = 1:n
    pii2sli[sli2pii[si]] = si
    pii2swi[swi2pii[si]] = si
  end
  SortedLinkedLW(l, w, sl, sw, sli2pii, swi2pii, pii2sli, pii2swi)
end

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

# TODO: check the possibility to use SortedLinkedLW to avoid creating
# tw new vectors every time
function becker2019_discretize(
  d :: Vector{D}, l :: Vector{S}, w :: Vector{S}, L :: S, W :: S;
  mark_single_piece = true
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
  mark_single_piece && (marks[l[1]] = true)
  y = l[1] # y: used to iterate capacity, inherited from knapsack papers
  for _ = 2:min(d[1], L ÷ l[1])
    marks[y += l[1]] = true
  end
  # If the lengths of the single pieces are not marked, the pairs are needed.
  if !mark_single_piece
    for pii = 2:N, pij = 1:(pii - 1) # PIece `j` (as opposed to `i`)
      mark[l[pii] + l[pij]] = true
    end
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
    # Mark cuts using just multiple copies of the same piece type (could be
    # done using a dummy cut at index zero with value true, but the array has
    # no index zero). 
    y = li
    mark_single_piece && (marks[y] = true)
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

function should_extract_piece_from_plate(
  pii :: D, L :: S, W :: S, symm :: UInt8, sllw :: SortedLinkedLW{D, S}
) :: Bool where {D, S}
  li = sllw.l[pii]
  wi = sllw.w[pii]
  (li > L || wi > W) && return false

  n = length(sllw.l)

  if symm == 1 || symm == 3
    fl = L - li
    fw = W

    for sli = 1:n
      sllw.sl[sli] > fl && break
      sllw.w[sllw.sli2pii[sli]] <= fw && return false
    end
  end

  if symm == 2 || symm == 3
    ll = L
    lw = W - wi

    for swi = 1:n
      sllw.sw[swi] > lw && break
      sllw.l[sllw.swi2pii[swi]] <= ll && return false
    end
  end

  return true
end

function should_extract_piece_from_plate(
  pii :: D, L :: S, W :: S, sllw :: SortedLinkedLW{D, S}
) :: Bool where {D, S}
  li = sllw.l[pii]
  wi = sllw.w[pii]
  (li > L || wi > W) && return false

  fl = L - li
  fw = W

  n = length(sllw.l)
  for sli = 1:n
    sllw.sl[sli] > fl && break
    sllw.w[sllw.sli2pii[sli]] <= fw && return false
  end

  ll = L
  lw = W - wi

  for swi = 1:n
    sllw.sw[swi] > lw && break
    sllw.l[sllw.swi2pii[swi]] <= ll && return false
  end

  return true
end

function partitions_no_symm(
  ::Type{P}, d :: Vector{D}, sllw :: SortedLinkedLW{D, S}, L :: S, W :: S
) where {D, S, P}
  l = sllw.l
  w = sllw.w
  @assert length(d) == length(l)
  @assert length(d) == length(w)
  l_not_single = becker2019_discretize(
    d, l, w, L, W; mark_single_piece = false
  ) where {D, S}
  max_piece_type = convert(D, length(l))
  max_num_plates = convert(P, L) * convert(P, W) * 3
  #dl = becker2019_discretize(d, l, L ÷ 2)
  #dw = becker2019_discretize(d, w, W ÷ 2)
  min_pil = minimum(l)
  min_piw = minimum(w)
  # If (n1, n2, n3) is in hcut (vcuts), then n1 may be partitioned in n2 and n3
  # by the means of a horizontal (vertical) cut. Not every possible partition
  # is present, just the ones which each child plate can fit at least one
  # piece, and the ones where one child plate has the same dimension as a piece
  # and the other is the dummy waste plate.
  hcuts = Vector{NTuple{3, P}}()
  vcuts = Vector{NTuple{3, P}}()
  # If the plate has the exact same size that a piece, then it is added to
  # the vector at pii2plis[piece_type].
  pii2plis = [Vector{P}() for _ = 1:max_piece_type]
  # If the piece has the exact same size that a plate, then it is added to
  # the vector at pli2piis[plate_type].
  pli2piis = Vector{D}[]
  # The list of plates attributes: plate length, plate width, and plate bound.
  # The plate index is the same as the index in pli_lwb.
  pli2lwsb = Vector{Tuple{S, S, UInt8, P}}()
  # plis: matrix of the plate dimensions in which zero means "never seen that
  # plate before" and nonzero means "this nonzero number is the plate index".
  plis = zeros(P, L, W, 3)
  plis[L, W, 3] = one(P)
  # next: plates already indexed but not yet processed, starts with (L, W,
  # 3, 1). The size of the original plate, the code for being able to cut
  # in any direction (same as third dimension of plis), and the plate index.
  # Storing the index as the third value is not necessary (as it could be
  # queried from plis) but this is probably more efficient this way.
  next = Vector{Tuple{S, S, UInt8, P}}()
  next_idx = one(P)
  #sizehint!(next, max_num_plates)
  push!(next, (L, W, 3, one(P)))
  # n: The amount of plates (the index of the highest plate type).
  n = one(P) # there is already the original plate
  while next_idx <= length(next)
    # PLate Length, Width, Symmetry, and Index
    pll, plw, pls, pli = next[next_idx]
    next_idx += 1
    plb = (L ÷ pll) * (W ÷ plw) # PLate Bound
    # It is not necessary to store the plate id in pli2lwsb because they are
    # added in order (a queue is used), so the array index is the plate index.
    push!(pli2lwsb, (pll, plw, pls, plb))
    # the pli2piis vector grows with the number of plates processed
    push!(pli2piis, D[])
    for pii in 1:max_piece_type # pii: PIece Index
      pil, piw = l[pii], w[pii] # PIece Length and Width (homophones, I know)
      if pil == pll && piw == plw
        # If the current plate is exactly the size of a plate, add it to
        # the list of plate codes that correspond to a piece.
        push!(pii2plis[pii], pli)
        push!(pli2piis[pli], pii)
      elseif should_extract_piece_from_plate(pii, pll, plw, pls, sllw)
        # If there is no way to cut the current piece from the current plate
        # except if we cut after the half of the plate... (i.e., the most
        # unbalanced cut possible, that is, making a cut with the minimal
        # length or width among the pieces, removes the possibility of getting
        # such piece from the larger child plate)
        if pil < pll && (pls == 1 || pls == 3)
          # If the current plate is allowed to be cut horizontally (i.e.,
          # reducing the length) and we need to reduce the length to make the
          # plate into the piece, then we do so.
          if iszero(plis[pil, plw, 2])
            push!(next, (pil, plw, 2, n += 1))
            plis[pil, plw, 2] = n
          end
          push!(hcuts, (pli, plis[pil, plw, 2], 0))
        elseif piw < plw && (pls == 2 || pls == 3)
          # If the current plate is allowed to be cut vertically (i.e.,
          # reducing the width) and we need to reduce the width to make the
          # plate into the piece, then we do so.
          if iszero(plis[pll, piw, 1])
            push!(next, (pll, piw, 1, n += 1))
            plis[pll, piw, 1] = n
          end
          push!(vcuts, (pli, plis[pll, piw, 1], 0))
        else
          #error("pls value is $(pls) but should be in [1,3]")
          println("happened")
        end
      end
    end
    if pls == 1 || pls == 3
      for y in becker2019_discretize(d, l, w, pll ÷ 2, plw)
        #y > pll ÷ 2 && break
        @assert plw >= min_piw
        @assert y >= min_pil
        @assert pll - y >= min_pil
        if iszero(plis[y, plw, 2])
          push!(next, (y, plw, 2, n += 1))
          plis[y, plw, 2] = n
        end
        if iszero(plis[pll - y, plw, 3])
          push!(next, (pll - y, plw, 3, n += 1))
          plis[pll - y, plw, 3] = n
        end
        push!(hcuts, (pli, plis[y, plw, 2], plis[pll - y, plw, 3]))
      end
    end
    if pls == 2 || pls == 3
      for x in becker2019_discretize(d, w, l, plw ÷ 2, pll)
        #x > plw ÷ 2 && break
        @assert pll >= min_pil
        @assert x >= min_piw
        @assert plw - x >= min_piw
        if iszero(plis[pll, x, 1])
          push!(next, (pll, x, 1, n += 1))
          plis[pll, x, 1] = n
        end
        if iszero(plis[pll, plw - x, 3])
          push!(next, (pll, plw - x, 3, n += 1))
          plis[pll, plw - x, 3] = n
        end
        push!(vcuts, (pli, plis[pll, x, 1], plis[pll, plw - x, 3]))
      end
    end
  end

  return pli2lwsb, hcuts, vcuts, pii2plis, pli2piis
end

function partitions(
  ::Type{P}, d :: Vector{D}, sllw :: SortedLinkedLW{D, S}, L :: S, W :: S
) where {D, S, P}
  l = sllw.l
  w = sllw.w
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
  # The list of plates attributes: plate length, plate width, and plate bound.
  # The plate index is the same as the index in pli_lwb.
  pli_lwb = Vector{Tuple{S, S, P}}()
  # plis: matrix of the plate dimensions in which zero means "never seen that
    # plate before" and nonzero means "this nonzero number is the plate index".
  plis = zeros(P, L, W)
  plis[L, W] = one(P)
  # next: plates already indexed but not yet processed, starts with (L, W, 1).
  # Storing the index as the third value is not necessary (as it could be
  # queried from plis) but this is probably more efficient this way.
  next = Vector{Tuple{S, S, P}}()
  next_idx = one(P)
  #sizehint!(next, max_num_plates)
  push!(next, (L, W, one(P)))
  # n: The amount of plates (the index of the highest plate type).
  n = one(P) # there is already the original plate
  while next_idx <= length(next)
    pll, plw, pli = next[next_idx] # PLate Length, Width, and Index
    next_idx += 1
    plb = (L ÷ pll) * (W ÷ plw) # PLate Bound
    # It is not necessary to store the plate id in pli_lwb because they are added
    # in order (a queue is used), so the array index is the plate index.
    push!(pli_lwb, (pll, plw, plb))
    for pii in 1:max_piece_type # pii: PIece Index
      if should_extract_piece_from_plate(pii, pll, plw, sllw)
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

  return pli_lwb, hnnn, vnnn, np
end

end # module

