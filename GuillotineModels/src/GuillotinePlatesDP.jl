module GuillotinePlatesDP

export SortedLinkedLW
export becker2019_discretize, gen_cuts, gen_cuts_sb

# Structure for keeping the piece lengths and widths both in the original order
# and sorted order, and allow to access the piece index or the width (length)
# counterpart when iterating the sorted lengths (widths).
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

# IMPORTANT: this constructor does not copy l and w, so any changes to l or w
# will reflect in it, silently and completely invalidating the structure.
# The first parameter is the type of integer to be used in the indexes.
# This constructor basically creates a sorted copy of both vectors, then
# create the four inverse indexes, and return them all in a single structure.
# This constructor is not implemente in the most efficient way but it is
# called one single time for each instance, so, irrelevant.
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

# TODO: break this method into two (like it was done for the flow model)
#       one for reducing the input, other to discretize only over one dimension
# Intelligent discretization method that: (1) takes demand into account;
# (2) takes the dimension that is not being discretized into account. It can
# also be used to discover which lengths are only shared by single pieces
# and not linear combinations of two or more pieces.
# d: piece demand (d)
# l: piece dimensions for dim that is being discretized (may be `l` or `w`)
# w: piece dimensions for dim that is NOT being discretized (may be `l` or `w`)
# L: plate dimension for dim that is being discretized (may be `L` or `W)`
# W: plate dimension for dim that is NOT being discretized (may be `L` or `W`)
# mark_single_piece: instead of discretizing the plate dimension, give all
#   the piece dimensions (for the dimension being discretized) that are shared
#   only by single pieces and not by two or more pieces (even of the same type)
#   combinations.
function becker2019_discretize(
  d :: Vector{D}, l :: Vector{S}, w :: Vector{S}, L :: S, W :: S;
  only_single_pieces = false, ignore_W = false, ignore_d = false
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
    if l[i] <= L && (ignore_W || w[i] <= W)
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
  !only_single_pieces && (marks[l[1]] = true)
  y = l[1] # y: used to iterate capacity, inherited from knapsack papers
  for _ = 2:(ignore_d ? L ÷ l[1] : min(d[1], L ÷ l[1]))
    marks[y += l[1]] = true
  end
  # If the lengths of the single pieces are not marked, the pairs are needed.
  if only_single_pieces
    for pii = 2:N, pij = 1:(pii - 1) # PIece `j` (as opposed to `i`)
      lij = l[pii] + l[pij]
      lij <= L && (marks[lij] = true)
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
        for r = 2:(ignore_d ? L : di)
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
    !only_single_pieces && (marks[y] = true)
    for _ = 2:(ignore_d ? L ÷ li : min(di, L ÷ li))
      marks[y += li] = true
    end
  end

  cuts = Vector{S}()
  sizehint!(cuts, L)
  if only_single_pieces
    # If we just want to check if the piece lengths cannot be obtained from
    # linear combinations of other pieces, we just iterate the unique piece
    # lengths in increasing order, and if they are not marked, then there is
    # no linear combination of other pieces that give the same length.
    usl = unique!(sort(l))
    for y in usl
      !marks[y] && push!(cuts, y)
    end
  else
    for (position, is_marked) in enumerate(marks)
      is_marked && push!(cuts, position)
    end
  end
  cuts
end

# Check if such piece should be extracted from a plate. This method may be used
# from both symmetry-breaking and non-symmetry-breaking contexts, because
# the default value of parameter symm is the allow-any-cut-direction value.
# A piece should be extracted from a plate if:
# (1) The plate allows the piece size.
# (2) The plate will not be cut in a way that leaves a second child large
#     enough to allow the piece size (i.e., if the piece extraction may be
#     postponed to be done in a second child of the current plate, it will be
#     postponed).
function should_extract_piece_from_plate(
  pii :: D, L :: S, W :: S, sllw :: SortedLinkedLW{D, S}, symm :: UInt8 = 0x03
) :: Bool where {D, S}
  li = sllw.l[pii]
  wi = sllw.w[pii]
  # If the piece does not fit the plate, then it is not even possible to
  # extract the piece from the plate.
  (li > L || wi > W) && return false

  n = length(sllw.l)

  # If we cut the plate by the length of the piece (a horizontal cut), does
  # this leaves space for another piece? If it does then we postpone to extract
  # the piece from a child of the current plate. Note that we can only do that
  # if the plate is marked to allow horizontal cuts or any cuts.
  # (The 'f' in the variables is from 'frontal'.)
  if symm == 1 || symm == 3
    fl = L - li
    fw = W

    for sli = 1:n # for the pieces in increase-or-stay order of their length
      # if the current piece does not fit the length none after it will fit
      sllw.sl[sli] > fl && break
      # if a piece fits the frontal child, then postpone extraction
      sllw.w[sllw.sli2pii[sli]] <= fw && return false
    end
  end

  if symm == 2 || symm == 3
    # The same as above but for the width of the piece (a vertical cut).
    # (The 'l' in the variables is from 'lateral'.)
    ll = L
    lw = W - wi

    # the same comments of the last loop are valid for this one
    for swi = 1:n
      sllw.sw[swi] > lw && break
      sllw.l[sllw.swi2pii[swi]] <= ll && return false
    end
  end

  # If no piece is small enough to make the plate-parameter be cut in a length
  # (or width) that leaves a second child large enough for extracting the
  # piece-parameter, then we cannot postpone and should allow extracting the
  # piece-parameter from the plate-parameter.
  return true
end

# TODO: internal change. The symmetry tag dimension of plis should be the
# first one, given how columnwise memory layout works.
function gen_cuts_sb(
  ::Type{P}, d :: Vector{D}, sllw :: SortedLinkedLW{D, S}, L :: S, W :: S;
  ignore_2th_dim = false, ignore_d = false
) where {D, S, P}
  (ignore_d || ignore_2th_dim) && @error "ignore_2th_dimm and ignore_d are not yet implemented for gen_cuts_sb, also, first improve the performance by memoizing the discretizations"
  l = sllw.l
  w = sllw.w
  @assert length(d) == length(l)
  @assert length(d) == length(w)
  only_single_l = becker2019_discretize(
    d, l, w, L, W; only_single_pieces = true, ignore_W = ignore_2th_dim,
    ignore_d = ignore_d
  )
  only_single_w = becker2019_discretize(
    d, w, l, W, L; only_single_pieces = true, ignore_W = ignore_2th_dim,
    ignore_d = ignore_d
  )
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
  piece_cuts_avoided_by_symmb = 0
  while next_idx <= length(next)
    # PLate Length, Width, Symmetry, and Index
    pll, plw, pls :: UInt8, pli = next[next_idx]
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
      elseif should_extract_piece_from_plate(pii, pll, plw, sllw, pls)
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
          piece_cuts_avoided_by_symmb += 1
        end
      end
    end
    # TODO: after is working, refactor the code below to use a single loop
    # using a vector saved in a variable, or an inner function
    sorted_in(a, v) = searchsortedfirst(a, v) <= length(a)
    if pls == 1 && sorted_in(only_single_w, plw)
      indexes = filter(i -> w[i] == plw && l[i] < pll, collect(1:max_piece_type))
      #@assert !isempty(indexes)
      for y in l[indexes]
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
    elseif pls == 1 || pls == 3
      for y in becker2019_discretize(
        d, l, w, pll ÷ 2, plw; ignore_2th_dim = ignore_2th_dim
      )
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
    if pls == 2 && sorted_in(only_single_l, pll)
      indexes = filter(i -> l[i] == pll && w[i] < plw, collect(1:max_piece_type))
      #@assert !isempty(indexes)
      for x in w[indexes]
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
    elseif pls == 2 || pls == 3
      for x in becker2019_discretize(
        d, w, l, plw ÷ 2, pll; ignore_2th_dim = ignore_2th_dim
      )
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
  #@show piece_cuts_avoided_by_symmb

  # For each possible subplatex size, same_size_plis groups the plate indexes
  # that refer to the plates with the same size, and distinct symmetry tags.
  same_size_plis = Vector{Vector{P}}()
  # The length of same_size_plis may be as low as one third of length(pli2lwsb)
  # and as high as length(pli2lwsb).
  sizehint!(same_size_plis, length(pli2lwsb))
  for j = 1:W, i = 1:L # column-major for performance
    temp = Vector{P}()
    !iszero(plis[i, j, 1]) && push!(temp, plis[i, j, 1])
    !iszero(plis[i, j, 2]) && push!(temp, plis[i, j, 2])
    !iszero(plis[i, j, 3]) && push!(temp, plis[i, j, 3])
    !iszero(length(temp)) && (push!(same_size_plis, temp); temp = Vector{P}())
  end
  @assert last(same_size_plis) == P[1]
  pop!(same_size_plis)

  return pli2lwsb, hcuts, vcuts, pii2plis, pli2piis, same_size_plis
end

function gen_cuts(
  ::Type{P}, d :: Vector{D}, sllw :: SortedLinkedLW{D, S}, L :: S, W :: S;
  ignore_2th_dim = false, ignore_d = false
) where {D, S, P}
  l = sllw.l
  w = sllw.w
  @assert length(d) == length(l)
  @assert length(d) == length(w)
  max_piece_type = convert(D, length(l))
  max_num_plates = convert(P, L) * convert(P, W)
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
  # Memoized discretizations. The discretized lengths for every plate width.
  # The discretized widths for every plate length.
  dls = [Vector{S}() for _ = 1:W]
  dws = [Vector{S}() for _ = 1:L]
  dl = dls[W] = becker2019_discretize(
    d, l, w, L, W; ignore_W = ignore_2th_dim, ignore_d = ignore_d
  )
  dw = dws[L] = becker2019_discretize(
    d, w, l, W, L; ignore_W = ignore_2th_dim, ignore_d = ignore_d
  )
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
    if !ignore_2th_dim
      isempty(dls[plw]) && (dls[plw] = becker2019_discretize(
        d, l, w, L, plw; ignore_d = ignore_d
      ))
      dl = dls[plw]
    end
    for y in dl
      y > pll ÷ 2 && break
      @assert plw >= min_piw
      @assert y >= min_pil
      @assert pll - y >= min_pil
      if iszero(plis[y, plw])
        push!(next, (y, plw, n += 1))
        plis[y, plw] = n
      end
      dl_sc_ix = searchsortedlast(dl, pll - y)
      @assert !iszero(dl_sc_ix)
      dl_sc = dl[dl_sc_ix]
      if iszero(plis[dl_sc, plw])
        push!(next, (dl_sc, plw, n += 1))
        plis[dl_sc, plw] = n
      end
      push!(hnnn, (pli, plis[y, plw], plis[dl_sc, plw]))
    end
    if !ignore_2th_dim
      isempty(dws[pll]) && (dws[pll] = becker2019_discretize(
        d, w, l, W, pll; ignore_d = ignore_d
      ))
      dw = dws[pll]
    end
    for x in dw
      x > plw ÷ 2 && break
      @assert pll >= min_pil
      @assert x >= min_piw
      @assert plw - x >= min_piw
      if iszero(plis[pll, x])
        push!(next, (pll, x, n += 1))
        plis[pll, x] = n
      end
      dw_sc_ix = searchsortedlast(dw, plw - x)
      @assert !iszero(dw_sc_ix)
      dw_sc = dw[dw_sc_ix]
      if iszero(plis[pll, dw_sc])
        push!(next, (pll, dw_sc, n += 1))
        plis[pll, dw_sc] = n
      end
      push!(vnnn, (pli, plis[pll, x], plis[pll, dw_sc]))
    end
  end

  return pli_lwb, hnnn, vnnn, np
end

end # module

