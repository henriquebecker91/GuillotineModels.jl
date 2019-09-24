module HeuristicsGC2D

# Returns zero if no slice 1:k (for any k) exists such that sum(p[1:k]) > bkv
# and sum(a[1:k]) <= A; otherwise, returns the smallest k value which satisfy
# these conditions.
function promising_kfirst(
  d :: Vector{D},
  p :: Vector{P},
  a :: Vector{P},
  A :: P,
  bkv :: P
) :: D where {P}
  n = length(p)
  @assert length(a) == n
  sump, suma = zero(P), zero(P)
  for k = 1:n
    sump += p[k]
    suma > A && return zero(D)
    suma += a[k]
    sump > bkv && return convert(D, k)
  end
  # the line below is only executed if sum(p) >= bkv
  return zero(D)
end

# Note: the rng object is modified, as it would be expected.
function iterated_greedy(
  d :: Vector{D}, p :: Vector{P}, l :: Vector{S}, w :: Vector{S},
  L :: S, W :: S, rng, max_iter_since_last_improv :: I = 1000000
) :: Tuple{P, Vector{D}} where {D, S, P, I}
  n = length(d) # number of distinct piece types
  A = convert(P, L) * convert(P, W)
  # assert explanation: all vectors should be the same length.
  @assert [n] == unique!(length.([p, l, w]))
  # loop control variables
  iter_last_improv = 0 # last iteration in which bkv changed
  curr_iter = 0 # current iteration
  bkv = 0 # best known value
  bks = zeros(D, n) # best known solution
  # Create a first ordering by efficiency (i.e., cost-benefit, profit-to-area
  # ratio). This way, the same instance with the same RNG will have the same
  # results, independent from the initial ordering of the piece types.
  a = @. convert(l, P) * convert(w, P)
  e = @. p / convert(Float64, a)
  o = sortperm(e)
  # remove any reference to the originals (the only object modified is the rng)
  d, p, l, w, a = deepcopy.((d, p, l, w, a))
  while curr_iter - iter_since_last_improv <= max_iter_since_last_improv
    curr_iter += 1
    permute!.((d, p, l, w, a), (o, o, o, o, o))
    k = promising_kfirst(d, p, a, A, bkv)
    if !iszero(k)
      @view d_, p_, l_, w_ = d[1:k], p[1:k], l[1:k], w[1:k]
      o_ = sortperm(w_; rev = true)
      permute!.((d_, p_, l_, w_), (o_, o_, o_, o_))
      v, s = first_fit_decr(d_, p_, l_, w_, L, W)
      if v > bkv
        bkv, bks = v, s
        iter_last_improv = curr_iter
      end
    end
    shuffle!(rng, o)
  end

  @assert sum(bks .* p) == bkv
  bkv, bks
end

# TODO: function needs to return solution structure to warm-start model
function first_fit_decr(
  d :: Vector{D}, p :: Vector{P}, l :: Vector{S}, w :: Vector{S},
  L :: S, W :: S
) :: Tuple{P, Vector{D}} where {D, S, P}
  n = length(d) # number of distinct piece types
  # assert explanation: all vectors should be the same length.
  @assert [n] == unique!(length.([p, l, w]))
  @assert issorted(w)
  lw = Vector{S}() # levels width
  ll = Vector{S}() # levels length
  rw = W # remaining width

  bkv = 0
  bks = zeros(D, n)
  for i in 1:n # for each piece type
    for _ in 1:d[i] # for each copy of the same piece type
      assigned = false
      for j in 1:length(ll) # for each already open level 
        # assert reasoning: the levels are created with the width of previous
        # pieces, and the pieces are ordered by decreasing weight, therefore
        # no piece may have more width than an already open level.
        @assert lw[j] >= w[i]
        if ll[j] >= l[i] # checks if the level has enough length yet
          ll[j] -= l[i]
          bkv += p[i]
          bks[i] += one(D)
          assigned = true
          break
        end
      end
      # If the piece does not fit any already open level...
      if !assigned
        # ... and there is no space to open a new level, then jump to the
        # next piece type.
        rw < w[i] && @goto next_piece_type
        # ... but there is yet space to open a new level, create a new level.
        push!(lw, w[i])
        push!(ll, L - l[i])
        bkv += p[i]
        bks[i] += one(D)
      end
    end
    @label next_piece_type
  end

  return bkv, bks
end

