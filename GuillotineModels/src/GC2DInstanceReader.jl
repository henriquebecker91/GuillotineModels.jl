module GC2DInstanceReader

export read_instance

# TODO: allow for out of order plate ids; save ids; 
# TODO: SOME INSTANCES ASSUME THE READER IS COMPLETELY INSENSTIVE TO
# WHITESPACE, CHANGE THE METHOD TO GET THE NEXT TOKEN AND NOT WORK
# ORIENTED BY LINES BUT BY WORDS/TOKENS/NUMBERS
function read_instance(filepath; config = Dict{Symbol,Bool}())
  default_config = Dict{Symbol,Bool}(
    :expect_pid => false
  )
  for k in keys(config)
    if (!haskey(default_config, k))
      @warn "the key " * string(k) * " is not recognized by " *
            "read_instance"
    end
  end
  # If the key appears in default_config and config, the merge keeps
  # the value in config.
  c = merge(default_config, config)
  (l, w, p, q) = (
    Vector{Int64}(), Vector{Int64}(), Vector{Int64}(), Vector{Int64}()
  )
  L = W = N = 0 # are set below but need to de defined here
  open(filepath) do f
    lnc = 1
    for ln in eachline(f)
      if      lnc == 1
        (L, W) = map(x->parse(Int64, x), split(ln))
      elseif  lnc == 2
        N = parse(Int64, ln)
      else    # all other lines
        # The instances of paper "Improved state space relaxation for
        # constrained two-dimensional guillotine cutting problems" have a N
        # different of the real number of lines because the instances with 50
        # and 25 items are the same, only changing the value (they could
        # have removed the unused lines...).
        #if (lnc > N + 2)
        #  @warn "In Olinto2019.read_instance, there are more plate-defining" *
        #    " lines than the given number of lines."
        #end
        # cp_ == current plate
        if c[:expect_pid]
          (cp_id, cp_l, cp_w, cp_p, cp_q) = map(x->parse(Int64, x), split(ln))
        else
          (cp_l, cp_w, cp_p, cp_q) = map(x->parse(Int64, x), split(ln))
        end
        if (c[:expect_pid] && cp_id != lnc - 2)
          @warn "In read_instance, line with id " * string(cp_id) *
            " but id expected was " * string(lnc - 2) * "."
        end
        push!(l, cp_l)
        push!(w, cp_w)
        push!(p, cp_p)
        push!(q, cp_q)
      end
      if lnc - 2 == N # use N and not all the lines in the file
        break
      end
      lnc = lnc + 1
    end
  end
  
  return L, W, l, w, p, q
end

end # module

