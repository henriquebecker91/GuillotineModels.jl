module GuillotinePlatesDP

using DelimitedFiles, StructArrays
push!(LOAD_PATH, "./")
using GC2DInstanceReader

export all_plates, nibnn, write_nibnn, partitions

function fits(l :: S, w :: S, L :: S, W :: S) :: Bool where {S}
  l <= L && w <= W
end

function the_four_children(
  l :: S, w :: S, L :: S, W :: S
) :: NTuple{4, NTuple{2, S}} where {S}
  # vertical frontal/lateral child
  vfc = (L - l, w)
  vlc = (L, W - w)
  # horizontal frontal/lateral child
  hfc = (L - l, W)
  hlc = (l, W - w)

  return (vfc, vlc, hfc, hlc)
end

function all_plates(
  l :: Vector{S}, w :: Vector{S}, L :: S, W :: S
) :: Vector{NTuple{2, S}} where {S}
  @assert length(l) == length(w)
  n = length(l)
  max_num_plates = L * W
  # seen: boolean matrix to check if a plate was already added to next
  seen = fill(false, L + 1, W + 1)
  # next: plates already seen but not yet processed, starts with (L, W)
  next = Vector{NTuple{2, S}}()
  sizehint!(next, max_num_plates)
  push!(next, (L, W))
  # used: plates that left 'next' and are able to fit at least one piece
  used = empty(next)
  sizehint!(used, max_num_plates)
  while !isempty(next)
    p_used = false
    pl, pw = pop!(next)
    for i in 1:n
      li, wi = l[i], w[i]
      fits(li, wi, pl, pw) ? p_used = true : continue
      # for every possible children plate
      for (cpl, cpw) in the_four_children(li, wi, pl, pw)
        seen[cpl + 1, cpw + 1] && continue
        seen[cpl + 1, cpw + 1] = true
        push!(next, (cpl, cpw))
      end
    end
    p_used && push!(used, (pl, pw))
  end

  return used
end

function nibnn(
  d :: Vector{D}, l :: Vector{S}, w :: Vector{S}, L :: S, W :: S
) :: Tuple{S, Vector{NTuple{5, S}}} where {D, S}
  @assert length(l) == length(w)
  max_piece_type = length(l)
  max_num_plates = L * W
  # ntag: matrix of the plate dimensions in which zero means "never seen that
  # plate before" and nonzero means "this nonzero number is the plate code".
  ntag = zeros(S, L + 1, W + 1)
  ntag[L + 1, W + 1] = one(S)
  # next: plates already tagged but not yet processed, starts with (L, W, 1).
  # Storing the tag as the third value is not necessary (as it could be queried
  # from ntag) but this is probably more efficient this way.
  next = Vector{NTuple{3, S}}()
  sizehint!(next, max_num_plates)
  push!(next, (L, W, one(S)))
  # used: plates that left 'next' and are able to fit at least one piece
  # (but their children can be waste) starts empty, as would be expected.
  # If a children is waste, its index (one of the last two values of the
  # 5-tuple) has a zero.
  used = Vector{NTuple{5, S}}()
  sizehint!(used, max_num_plates)
  n = one(S) # there is already the original plate
  while !isempty(next)
    pll, plw, plt = pop!(next) # PLate Length, Width, and Tag
    how_many_of_this_plate_fit_original = (L ÷ pll) * (W ÷ plw)
    p_used = false # If the plate has at least one piece that fits it.
    for pii in 1:max_piece_type # pii: PIece Index
      pil, piw = l[pii], w[pii] # PIece Length and Width (homophones, I know)
      # If this piece does not fit the current plate skip to the next.
      fits(pil, piw, pll, plw) ? p_used = true : continue
      # Vertical/Horizontal Frontal/Lateral Children
      vfc, vlc, hfc, hlc = the_four_children(pil, piw, pll, plw)
      # The combinations of `n`, `i`, and `b` are always unique because `n` is
      # unique (`ntag` guarantees each plate type `n` is only processed a
      # single time), however, this triple combination can waste one or both
      # children (zero value) in the two possible cuts, and two (n, i, b, 0, 0)
      # (or (n, i, b, 0, lct), (n, i, b, fct, 0)) will be added.
      # The variables below avoid that.
      second_time = false
      previous_ftag = previous_ltag = zero(S)
      # Note that the loop below will execute two times, one for both vertical
      # childs, and another for both horizontal childs. The code is the same
      # independent from which cut was done.
      # Frontal/Lateral Child Length/Width
      for ((fcl, fcw), (lcl, lcw)) in ((vfc, vlc), (hfc, hlc))
        ftag = ntag[fcl + 1, fcw + 1] # Frontal child TAG
        # If the child plate was not seen before, check if it should be tagged
        # (i.e., it is a plate that can fit at least one piece inside it).
        if iszero(ftag) && any(((a, b),) -> fits(a, b, fcl, fcw), zip(l, w))
          # Increments `n`, saves the plate as tagged, and add it to `next`.
          ftag = ntag[fcl + 1, fcw + 1] = (n += one(S))
          push!(next, (fcl, fcw, ftag))
        end
        ltag = ntag[lcl + 1, lcw + 1] # Lateral child TAG
        # See the comments for the first `if` of this loop.
        if iszero(ltag) && any(((a, b),) -> fits(a, b, lcl, lcw), zip(l, w))
          ltag = ntag[lcl + 1, lcw + 1] = (n += one(S))
          push!(next, (lcl, lcw, ltag))
        end
        # To understand the `if` below check the comments above the
        # definition of the variables used in the `if`.
        if second_time && ftag == previous_ftag && ltag == previous_ltag
          # Note that second_time guarantees this is the second loop
          # so break or continue are the same.
          break 
        end
        # Here we finally add a combination of the five indexes to used, this
        # is done by the min between "the number of times the plate can appear
        # on the original plate" and the "max number of times we would cut pii
        # from some plate" .
        for i = 1:min(how_many_of_this_plate_fit_original, d[pii])
          # plt = plate tag; i = i-esim plate of the type; pii = the piece
          # index cut from this plate; ftag/ltag the tag of the children
          # plates (can be zero if the plate is waste).
          push!(used, (plt, i, pii, ftag, ltag))
        end
        second_time = true
        previous_ltag = ltag
        previous_ftag = ftag
      end
    end
    # NOTE: except by the number one plate, every other plate is only added
    # to `next` if it can fit at least one piece inside it, so p_used has to
    # be marked as true at the end of the loop, or we did something wrong.
    # And yes, it is only used for this assert, it is not used anywhere else.
    @assert p_used
  end

  return n, used
end

function write_nibnn(src, dest)
  L, W, l, w, _, d = read_instance(src)
  open(dest, "w") do io
    n, nibnn_tuples = nibnn(d, l, w, L, W)
    println(io, string(n))
    nibnn_matrix = hcat(fieldarrays(StructArray(nibnn_tuples))...)
    writedlm(io, nibnn_matrix, ' ')
  end
end

function write_nibnn(
  d :: Vector{D}, l :: Vector{S}, w :: Vector{S}, L :: S, W :: S, fname 
) where {D, S}
  open(fname, "w") do io
    n, nibnn_tuples = nibnn(d, l, w, L, W)
    println(io, string(n))
    nibnn_matrix = hcat(fieldarrays(StructArray(nibnn_tuples))...)
    writedlm(io, nibnn_matrix, ' ')
  end
end

function partitions(
  ::Type{D}, ::Type{P}, l :: Vector{S}, w :: Vector{S}, L :: S, W :: S
) :: Tuple{P, Vector{Tuple{P, D, P, P}}} where {D, S, P}
  @assert length(l) == length(w)
  max_piece_type = convert(D, length(l))
  max_num_plates = convert(P, L) * convert(P, W)
  # plis: matrix of the plate dimensions in which zero means "never seen that
  # plate before" and nonzero means "this nonzero number is the plate index".
  plis = zeros(P, L + 1, W + 1)
  plis[L + 1, W + 1] = one(P)
  # next: plates already indexed but not yet processed, starts with (L, W, 1).
  # Storing the index as the third value is not necessary (as it could be
  # queried from plis) but this is probably more efficient this way.
  next = Vector{Tuple{S, S, P}}()
  #sizehint!(next, max_num_plates)
  push!(next, (L, W, one(P)))
  # used: plates that left 'next' and are able to fit at least one piece
  # (but their children can be waste). Starts empty, as would be expected.
  # If a child is waste, its index (one of the last two values of the
  # triple) is zero.
  used = Vector{Tuple{P, D, P, P}}()
  #sizehint!(used, max_num_plates)
  n = one(P) # there is already the original plate
  while !isempty(next)
    pll, plw, pli = pop!(next) # PLate Length, Width, and Tag
    p_used = false # If the plate has at least one piece that fits it.
    for pii in 1:max_piece_type # pii: PIece Index
      pil, piw = l[pii], w[pii] # PIece Length and Width (homophones, I know)
      # If this piece does not fit the current plate skip to the next.
      fits(pil, piw, pll, plw) ? p_used = true : continue
      # Vertical/Horizontal Frontal/Lateral Children
      vfc, vlc, hfc, hlc = the_four_children(pil, piw, pll, plw)
      # The combinations of `n` and `b` are always unique because `n` is
      # unique (`plis` guarantees each plate type `n` is only processed a
      # single time), however, this triple combination can waste one or both
      # children (zero value) in the two possible cuts, and two (n, i, b, 0, 0)
      # (or (n, i, b, 0, lct), (n, i, b, fct, 0)) will be added.
      # The variables below avoid that.
      second_time = false
      previous_fci = previous_lci = zero(S)
      # Note that the loop below will execute two times, one for both vertical
      # childs, and another for both horizontal childs. The code is the same
      # independent from which cut was done.
      # Frontal/Lateral Child Length/Width
      for ((fcl, fcw), (lcl, lcw)) in ((vfc, vlc), (hfc, hlc))
        fci = plis[fcl + 1, fcw + 1] # Frontal Child Index
        # If the child plate was not seen before, check if it should be
        # indexed (i.e., it is a plate that can fit at least one piece inside
        # it).
        if iszero(fci) && any(((a, b),) -> fits(a, b, fcl, fcw), zip(l, w))
          # Increments `n`, id the plate, add it to `next`, and increase
          # `used` size.
          fci = plis[fcl + 1, fcw + 1] = (n += one(S))
          push!(next, (fcl, fcw, fci))
        end
        lci = plis[lcl + 1, lcw + 1] # Lateral Child Index
        # See the comments for the first `if` of this loop.
        if iszero(lci) && any(((a, b),) -> fits(a, b, lcl, lcw), zip(l, w))
          lci = plis[lcl + 1, lcw + 1] = (n += one(S))
          push!(next, (lcl, lcw, lci))
        end
        # To understand the `if` below check the comments above the
        # definition of the variables used in the `if`.
        if second_time && fci == previous_fci && lci == previous_lci
          # Note that second_time guarantees this is the second loop
          # so break or continue are the same.
          break 
        end
        # pli = plate id; i = i-esim plate of the type; pii = the piece
        # index cut from this plate; fci/lci the id of the children
        # plates (can be zero if the plate is waste).
        push!(used, (pli, pii, fci, lci))

        second_time = true
        previous_lci = lci
        previous_fci = fci
      end
    end
    # NOTE: except by the number one plate, every other plate is only added
    # to `next` if it can fit at least one piece inside it, so p_used has to
    # be marked as true at the end of the loop, or we did something wrong.
    # And yes, it is only used for this assert, it is not used anywhere else.
    @assert p_used
  end

  return n, used
end

end # module

