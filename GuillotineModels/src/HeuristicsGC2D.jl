module HeuristicsGC2D

using Random

export first_fit_decr, iterated_greedy

# Returns zero if no slice 1:k (for any k, but without changing the ordering)
# exists such that sum(p[1:k]) > bkv and sum(a[1:k]) <= A; otherwise, returns
# the smallest k value which satisfy these conditions.
function promising_kfirst(
  d :: AbstractVector{D},
  p :: AbstractVector{P},
  a :: AbstractVector{P},
  A :: P,
  bkv :: P
) :: D where {D, P}
  n = length(d)
  @assert [n] == unique!(length.([p, a]))
  sump, suma = zero(P), zero(P)
  for k = 1:n
    suma += a[k]
    suma > A && return zero(D)
    sump += p[k]
    sump > bkv && return convert(D, k)
  end
  # the line below is only executed if sum(p) <= bkv
  return zero(D)
end

function expand(
  d :: AbstractVector{D}, a :: AbstractVector{T}
) :: AbstractVector{T} where {D, T}
  n = length(d)
  @assert length(a) == n
  collect(Base.Iterators.Flatten(
    map(
      ((v, q),) -> repeat([v], q),
      zip(a, d)
    )
  ))
end

# Note: the rng object is modified, as it would be expected.
function iterated_greedy(
  d :: AbstractVector{D}, p :: AbstractVector{P}, l :: AbstractVector{S}, w :: AbstractVector{S},
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
  a = @. convert(P, l) * convert(P, w)
  e = @. p / convert(Float64, a)
  o = sortperm(e)
  i = collect(1:n) # permutation vector used to keep the original order
  # TODO: using i is not the most efficient way of dealing with the problem

  op = p # save original p for assert at end

  # remove any reference to the originals (the only object modified is the rng)
  d, p, l, w, a = deepcopy.((d, p, l, w, a))
  while curr_iter - iter_last_improv <= max_iter_since_last_improv
    curr_iter += 1
    permute!.((d, p, l, w, a, i), (o, o, o, o, o, o))
    k = promising_kfirst(d, p, a, A, bkv)
    if !iszero(k)
      @views d_, p_, l_, w_, a_, i_ = d[1:k], p[1:k], l[1:k], w[1:k], a[1:k], i[1:k]
      o_ = sortperm(w_; rev = true)
      permute!.((d_, p_, l_, w_, a_, i_), (o_, o_, o_, o_, o_, o_))
      v, s = first_fit_decr(d_, p_, l_, w_, L, W)
      if v > bkv
        (s_ = zeros(D, n))[1:k] = s
        invpermute!(s_, i)
        bkv, bks = v, s_
        iter_last_improv = curr_iter
      end
    end
    shuffle!(rng, o)
  end

  @assert sum(bks .* op) == bkv
  bkv, bks
end

# Version that sorts the arrays by width before starting.
function first_fit_decr!(
  d :: AbstractVector{D}, p :: AbstractVector{P}, l :: AbstractVector{S}, w :: AbstractVector{S},
  L :: S, W :: S
) :: Tuple{P, AbstractVector{D}} where {D, S, P}
  o = sortperm(w; rev = true)
  permute!.((d, p, l, w), (o, o, o, o))
  return first_fit_decr(d, p, l, w, L, W)
end

# TODO: function needs to return solution structure to warm-start model
function first_fit_decr(
  d :: AbstractVector{D}, p :: AbstractVector{P}, l :: AbstractVector{S}, w :: AbstractVector{S},
  L :: S, W :: S
) :: Tuple{P, Vector{D}} where {D, S, P}
  n = length(d) # number of distinct piece types
  # assert explanation: all vectors should be the same length.
  @assert [n] == unique!(length.([p, l, w]))
  @assert issorted(w; rev = true)
  lw = Vector{S}() # levels width
  ll = Vector{S}() # levels length
  rw = W # remaining width

  bkv = 0
  bks = zeros(D, n)
  for i in 1:n # for each piece type
    for _ in 1:d[i] # for each copy of the same piece type
      assigned = false
      for j in 1:length(lw) # for each already open level 
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
        rw -= w[i]
        bkv += p[i]
        bks[i] += one(D)
      end
    end
    @label next_piece_type
  end

  @assert sum(bks .* p) == bkv
  return bkv, bks
end

end # module

