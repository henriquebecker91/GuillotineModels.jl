# Check if the necessary packages are installed and define the project
# to make them available to importation by "using".
using Pkg
Pkg.activate("..")
Pkg.instantiate()

# Load the necessary packages. The ones after the `push!` are from this
# repository.
using JuMP, ArgParse, MathOptInterface
#using CPLEX
#using_cplex = true
using_cplex = false
#using Cbc
using Gurobi
#using Profile

using MathOptFormat # disabled because dependency hell
using Random # for MersenneTwister used in heuristic
const MOI = MathOptInterface
push!(LOAD_PATH, "../src/")
using AllSubplatesModel, FlowModel, GC2DInstanceReader, HeuristicsGC2D
using ModelUtils # utility methods for dealing with JuMP models

# Create a new model with a configured underlying solver.
function new_empty_and_configured_model(p_args; no_solver_out = no_solver_out)
  @assert isone(length(p_args["threads"]))
  threads = first(p_args["threads"])
  @assert isone(length(p_args["seed"]))
  seed = first(p_args["seed"])
  @assert isone(length(p_args["det-time-limit"]))
  det_time_limit = first(p_args["det-time-limit"])
  @assert isone(length(p_args["time-limit"]))
  time_limit = first(p_args["time-limit"])

  if p_args["disable-solver-tricks"]
    m = JuMP.direct_model(
      #=CPLEX.Optimizer(
        CPX_PARAM_EPGAP = 1e-6,
        CPX_PARAM_PROBE = -1,
        CPX_PARAM_HEURFREQ = -1,
        CPX_PARAM_REPEATPRESOLVE = -1,
        CPX_PARAM_DIVETYPE = 1,
        CPX_PARAM_DETTILIM = det_time_limit,
	CPX_PARAM_TILIM = time_limit,
        #CPX_PARAM_VARSEL = CPLEX.CPX_VARSEL_MAXINFEAS,
        CPX_PARAM_STARTALG = CPLEX.CPX_ALG_BARRIER,
        CPX_PARAM_SCRIND = no_solver_out ? CPLEX.CPX_OFF : CPLEX.CPX_ON,
        CPX_PARAM_THREADS = threads,
        CPX_PARAM_RANDOMSEED = seed
      )=#
      Gurobi.Optimizer(
        Method = 2, # use barrier for LP
        #PreSparsify = 1, # try to reduce nonzeros
        #Presolve = 2, # aggressive presolving
        Threads = threads,
        Seed = seed,
        OutputFlag = no_solver_out ? 0 : 1,
        MIPGap = 1e-6,
        TimeLimit = time_limit
      )
      #=Cbc.Optimizer(
        threads = 1,
        ratioGap = 1e-6,
        logLevel = no_solver_out ? 0 : 1,
        randomSeed = seed,
        barrier = true,
	seconds = time_limit
      )=#
    )
  else
    m = JuMP.direct_model(
      #=CPLEX.Optimizer(
        CPX_PARAM_EPGAP = 1e-6,
        CPX_PARAM_DETTILIM = det_time_limit,
	CPX_PARAM_TILIM = time_limit,
        #CPX_PARAM_VARSEL = CPLEX.CPX_VARSEL_MAXINFEAS,
        CPX_PARAM_STARTALG = CPLEX.CPX_ALG_BARRIER,
        CPX_PARAM_SCRIND = no_solver_out ? CPLEX.CPX_OFF : CPLEX.CPX_ON,
        CPX_PARAM_THREADS = threads,
        CPX_PARAM_RANDOMSEED = seed
      )=#
      Gurobi.Optimizer(
        Method = 2, # use barrier for LP
        #PreSparsify = 1, # try to reduce nonzeros
        #Presolve = 2, # aggressive presolving
        Threads = threads,
        Seed = seed,
        OutputFlag = no_solver_out ? 0 : 1,
        MIPGap = 1e-6,
        TimeLimit = time_limit
      )
      #=Cbc.Optimizer(
        threads = 1,
        ratioGap = 1e-6,
        logLevel = no_solver_out ? 0 : 1,
        randomSeed = seed,
        barrier = true,
	seconds = time_limit
      )=#
    )
  end

  return m
end

# Read the instance, build the model, solve the model, and print related stats.
function read_build_solve_and_print(
  p_args, instfname_idx; no_csv_out = false, no_solver_out = false
)
  instfname = p_args["instfnames"][instfname_idx]
  if !no_csv_out
    @show instfname
    seed = first(p_args["seed"])
    @show seed
  end
  L, W, l, w, p, d = GC2DInstanceReader.read_instance(instfname)
  if p_args["div-and-round-nearest"][1] != 1
    factor = p_args["div-and-round-nearest"][1]
    roundmode = RoundNearest
  elseif p_args["div-and-round-up"][1] != 1
    factor = p_args["div-and-round-up"][1]
    roundmode = RoundUp
  elseif p_args["div-and-round-down"][1] != 1
    roundmode = RoundDown
    factor = p_args["div-and-round-down"][1]
  else
    factor = 1
  end
  if factor != 1
    S = eltype(L)
    L = convert(S, round(L / factor, roundmode))
    W = convert(S, round(W / factor, roundmode))
    l = convert.(S, round.(l ./ factor, roundmode))
    w = convert.(S, round.(w ./ factor, roundmode))
  end
  before_model_build = time()
  m = new_empty_and_configured_model(p_args; no_solver_out = no_solver_out)
  p_args["flow-model"] && p_args["only-binary-variables"] && @warn "the algorithm flow-model does not support the option only-binary-variables; ignoring the flag only-binary-variables"
  p_args["flow-model"] && p_args["break-hvcut-symmetry"] && @warn "the algorithm flow-model does not support the option break-hvcut-symmetry; ignoring the flag break-hvcut-symmetry"
  if p_args["flow-model"]
    FlowModel.build_model(m, d, p, l, w, L, W)
  else
    if p_args["break-hvcut-symmetry"]
      _, hvcuts, pli2lwsb, _, _ = AllSubplatesModel.build_model_with_symmbreak(
        m, d, p, l, w, L, W;
        only_binary = p_args["only-binary-variables"],
        use_c25 = p_args["use-c25"],
        ignore_2th_dim = p_args["ignore-2th-dim"],
        ignore_d = p_args["ignore-d"],
        round2disc = p_args["round2disc"]
      )
    else
      _, hvcuts, pli_lwb, np = AllSubplatesModel.build_model_no_symmbreak(
        m, d, p, l, w, L, W;
        only_binary = p_args["only-binary-variables"],
        use_c25 = p_args["use-c25"],
        ignore_2th_dim = p_args["ignore-2th-dim"],
        ignore_d = p_args["ignore-d"],
        round2disc = p_args["round2disc"],
        no_cut_position = p_args["no-cut-position"],
        no_redundant_cut = p_args["no-redundant-cut"],
        no_furini_symmbreak = p_args["no-furini-symmbreak"],
        faithful2furini2016 = p_args["faithful2furini2016"],
        lb = get(p_args["lower-bounds"], instfname_idx, 0),
        ub = get(p_args["upper-bounds"], instfname_idx, sum(d .* p)),
	print_debug = !no_csv_out
      )
    end
  end
  time_to_build_model = time() - before_model_build
  # The three lines below create a .mps file of the model before solving it.
  #=
  mps_model = MathOptFormat.MPS.Model()
  MOI.copy_to(mps_model, backend(m))
  MOI.write_to_file(mps_model, "last_run.mps")
  =#

  if !no_csv_out
    @show time_to_build_model
    n = length(d)
    @show n
    n_ = sum(d)
    @show n_
    p_ = max(L/minimum(l), W/minimum(w))
    @show p_
    num_vars = num_variables(m)
    @show num_vars
    num_constrs = num_all_constraints(m)
    @show num_constrs
    #@show m # just in case there is something here I did not output
    #print(m) # just in case there is something here I did not output
    if !p_args["flow-model"] && !p_args["break-hvcut-symmetry"]
      num_plates = length(pli_lwb)
      @show num_plates
      ps = m[:picuts]
      num_picuts = length(ps)
      @show num_picuts
      cm = m[:cuts_made]
      num_cuts_made = length(cm)
      @show num_cuts_made
    end
  end
  flush(stdout) # guarantee all messages will be flushed before calling cplex

  if p_args["warm-start"]
    @assert !p_args["faithful2furini2016"] && !p_args["flow-model"]
    (heur_obj, heur_sel, heur_pat), heur_time, _, _, _ = @timed iterated_greedy(
      d, p, l, w, L, W, MersenneTwister(p_args["seed"][1])
    )
    if !no_csv_out
      @show heur_time
      @show heur_obj
      @show heur_sel
      @show heur_pat
    end
    saved_bounds = AllSubplatesModel.disable_unrestricted_cuts!(
      m, sort(l), sort(w), hvcuts, pli_lwb
    )
    heur_ws_time = @elapsed AllSubplatesModel.warm_start(
      m, l, w, L, W, heur_pat, pli_lwb, hvcuts, np;
      faithful2furini2016 = p_args["faithful2furini2016"]
    )
    !no_csv_out && @show heur_ws_time
    restricted_time = @elapsed optimize!(m)
    !no_csv_out && @show restricted_time
    restricted_sol = value.(all_variables(m))
    restricted_objval = objective_value(m)
    AllSubplatesModel.restore_bound!.(saved_bounds)
    set_start_value.(all_variables(m), restricted_sol)
  end

  p_args["do-not-solve"] && return nothing
  #@error "testing furini, if optimize is called, all memory will be consumed"
  #exit(0)
  vars_before_deletes = all_variables(m)
  if p_args["final-pricing"] || p_args["relax2lp"]
    original_settings = relax_all_vars!(m)
  end
  time_to_solve_model = @elapsed optimize!(m)
  if (p_args["final-pricing"] || p_args["relax2lp"]) && !no_csv_out
    println("time_to_solve_relaxed_model = $(time_to_solve_model)")
  end
  if p_args["final-pricing"]
    # get the given lower bound on the instance if there is one
    best_lb = get(p_args["lower-bounds"], instfname_idx, 0)
    # if warm-start is enabled and the lb is better, use it
    p_args["warm-start"] && (best_lb = max(best_lb, restricted_objval))
    # delete all variables which would only reduce obj to below the lower bound
    #if !no_csv_out
    #  @profile (mask = delete_vars_by_pricing!(m, Float64(best_lb)))
    #  Profile.print()
    #else
    #  mask = delete_vars_by_pricing!(m, Float64(best_lb))
    #end
    kept = which_vars_to_keep(m, Float64(best_lb))
    unrelax_vars!(vars_before_deletes, m, original_settings)
    saved_bounds = fix_vars!(all_variables(m)[.!kept])
    if !no_csv_out
      num_vars_after_pricing = sum(kept)
      num_vars_del_by_pricing = length(kept) - num_vars_after_pricing
      @show num_vars_after_pricing
      @show num_vars_del_by_pricing
    end
    time_to_solve_model += @elapsed optimize!(m)
  end
  sleep(0.01) # shamefully needed to avoid out of order messages from cplex

  if !no_csv_out
    @show time_to_solve_model
    if primal_status(m) == MOI.FEASIBLE_POINT
      obj_value = objective_value(m)
    else
      obj_value = 0
    end
    @show obj_value
    if !p_args["relax2lp"]
      obj_bound = objective_bound(m)
      @show obj_bound
    end
    stop_reason = termination_status(m)
    @show stop_reason
    stop_code = Int64(stop_reason)
    @show stop_code
    if p_args["flow-model"]
      ps = value.(m[:edge][pii] for pii = 1:length(d) if is_valid(m, m[:edge][pii]))
      ps_nz = [iv for iv in enumerate(ps) if iv[2] > 0.001]
      @show ps_nz
    else
      if p_args["break-hvcut-symmetry"]
        ps = m[:pieces_sold]
        cm = m[:cuts_made]
        ps_nz = [(i, value(ps[i])) for i = 1:length(ps) if is_valid(m, ps[i]) && value(ps[i]) > 0.001]
        cm_nz = [(i, value(cm[i])) for i = 1:length(cm) if is_valid(m, cm[i]) && value(cm[i]) > 0.001]
        @show ps_nz
        @show cm_nz
        println("(piece length, piece width) => (piece index, amount in solution, profit of single piece, total profit contributed in solution) ")
        foreach(ps_nz) do e
          pii, v = e
          println((l[pii], w[pii]) => (pii, v, p[pii], v * p[pii]))
        end
        println("(parent plate length, parent plate width) => ((first child plate length, first child plate width), (second child plate length, second child plate width))")
        foreach(cm_nz) do e
          i, _ = e
          parent, fchild, schild = hvcuts[i]
          if iszero(schild)
            println((pli2lwsb[parent][1], pli2lwsb[parent][2]) => ((pli2lwsb[fchild][1], pli2lwsb[fchild][2]), (0, 0)))
          else
            println((pli2lwsb[parent][1], pli2lwsb[parent][2]) => ((pli2lwsb[fchild][1], pli2lwsb[fchild][2]), (pli2lwsb[schild][1], pli2lwsb[schild][2])))
          end
        end
      else # same as: if !p_args["break-hvcut-symmetry"]
        p_args["relax2lp"] && @warn "relax2lp flag used, be careful regarding solution"
        ps_nz = Vector{Tuple{Int, Float64}}()
        cm_nz = Vector{Tuple{Int, Float64}}()
        for i = 1:length(ps)
          if is_valid(m, ps[i]) && value(ps[i]) > 0.0001
            push!(ps_nz, (i, value(ps[i])))
          end
        end
        for i = 1:length(cm)
          if is_valid(m, cm[i]) && value(cm[i]) > 0.0001
            push!(cm_nz, (i, value(cm[i])))
          end
        end
        #ps_nz = [(i, value(ps[i])) for i = 1:length(ps) if (is_valid(m, ps[i]) && value(ps[i]) > 0.001)]
        #cm_nz = [(i, value(cm[i])) for i = 1:length(cm) if (is_valid(m, cm[i]) && value(cm[i]) > 0.001)]
        @show ps_nz
        @show cm_nz
        println("(plate length, plate width) => (number of times this extraction happened, piece length, piece width)")
        foreach(ps_nz) do e
          i, v = e
          pli, pii = np[i]
          println((pli_lwb[pli][1], pli_lwb[pli][2]) => (v, l[pii], w[pii]))
        end
        println("(parent plate length, parent plate width) => (number of times this cut happened, (first child plate length, first child plate width), (second child plate length, second child plate width))")
        foreach(cm_nz) do e
          i, v = e
          parent, fchild, schild = hvcuts[i]
	  @assert !iszero(parent)
	  @assert !iszero(fchild)
	  if iszero(schild)
	    @assert p_args["faithful2furini2016"]
	    if pli_lwb[parent][1] == pli_lwb[fchild][1]
              println((pli_lwb[parent][1], pli_lwb[parent][2]) => (v, (pli_lwb[fchild][1], pli_lwb[fchild][2]), (pli_lwb[parent][1], pli_lwb[parent][2] - pli_lwb[fchild][2])))
	    else
	      @assert pli_lwb[parent][2] == pli_lwb[fchild][2]
              println((pli_lwb[parent][1], pli_lwb[parent][2]) => (v, (pli_lwb[fchild][1], pli_lwb[fchild][2]), (pli_lwb[parent][1] - pli_lwb[fchild][1], pli_lwb[parent][2])))
	    end
	  else
            println((pli_lwb[parent][1], pli_lwb[parent][2]) => (v, (pli_lwb[fchild][1], pli_lwb[fchild][2]), (pli_lwb[schild][1], pli_lwb[schild][2])))
	  end
        end
      end
    end
  end

  return nothing
end

# Definition of the command line arguments.
function parse_script_args(args = ARGS)
    s = ArgParseSettings()
    @add_arg_table s begin
      "--threads"
        help = "number of threads used by the solver (not model building)"
        arg_type = Int
        default = [1]
        nargs = 1
      "--time-limit"
        help = "time limit in seconds for solver B&B (not model building, not root solving), default one year"
        arg_type = Float64
        default = [31536000.0]
        nargs = 1
      "--det-time-limit"
        help = "deterministic time limit for CPLEX B&B (not model building, not root solving), only CPLEX has deterministic time limit among the solvers"
        arg_type = Float64
        default = [1.0E+75]
        nargs = 1
      "--seed"
        help = "seed used by the solver (model building is deterministic)"
        arg_type = Int
        default = [1]
        nargs = 1
      "--disable-solver-tricks"
        help = "disable heuristics, probe, repeat presolve and dive"
        nargs = 0
      "--disable-solver-output"
        help = "disable the CPLEX output for all instances"
        nargs = 0
      "--do-not-solve"
        help = "just builds the model, do not call the solver over it"
        nargs = 0
      "--do-not-mock-first-for-compilation"
        help = "if passed then the first instance will not be solved twice"
        nargs = 0
      "--break-hvcut-symmetry"
        help = "will double the number of variables, and reduce symmetry"
        nargs = 0
      "--only-binary-variables"
        help = "CAUTION DO NOT GUARANTEE OPTIMALITY FOR NOW, IS HEURISTIC"
        nargs = 0
      "--use-c25"
        help = "add the tightening constraints 2.5 (ignored by flow)"
        nargs = 0
      "--round2disc"
        help = "round the second child size to a discretized position"
        nargs = 0
      "--ignore-2th-dim"
        help = "ignore the dimension not being discretized during discretization, used to measure impact (does not affect flow)"
        nargs = 0
      "--final-pricing"
        help = "uses the best lb available (from --lower-bounds or --warm-start) to remove variables after solving the continuous relaxation"
        nargs = 0
      "--faithful2furini2016"
        help = "tries to be the most faithful possible to the description on the Furini2016 paper (more specifically the complete model, with the reductions, FOR NOW without a heuristic solution first and without pricing); the flags --no-cut-position, --no-redundant-cut and --no-furini-symmbreak disable parts of this reimplementation"
        nargs = 0
      "--no-redundant-cut"
        help = "disables Furini2016 Redundant-Cut reduction: a bunch of flags is used to check if some trim cut is necessary for optimality, or is dominated by other trim cuts"
        nargs = 0
      "--no-furini-symmbreak"
        help = "disables Furini2016 symmetry breaking: as trim cuts are needed in faithful2furini2016 mode, the symmetry-breaking is less restrictive than 'do not cut after midplate', it just removes each cut y > midplate in which plate_size - y is an existing cut; enabling this flag makes the code just use all discretized positions as cuts"
        nargs = 0
      "--no-cut-position"
        help = "disables Furini2016 Cut-Position reduction: if a trivial heuristic shows that no combination of six or more pieces fit a plate, then the plate may be cut with restricted cuts without loss to optimality"
        nargs = 0
      "--warm-start"
        help = "(works only for Revised Furini) uses the heuristic described in Furini2016 to generate a initial primal feasible solution, warm-start the model with unrestricted cuts fixed to zero, and then unfix the unrestricted cuts to solve the complete model"
        nargs = 0
      "--relax2lp"
        help = "(works on Furini-like and Flow1) integer and binary variables become continuous"
        nargs = 0
      "--ignore-d"
        help = "ignore the demand information during discretization, used to measure impact (does not affect flow)"
        nargs = 0
      "--flow-model"
        help = "use the flow model instead of subplate model (ignore the --break-hvcut-symmetry and --only-binary-variables flags)"
        nargs = 0
      "--lower-bounds"
        help = "takes a single string, the string has N comma-separated-numbers where N is the number of instances passed"
        nargs = 1
      "--upper-bounds"
        help = "takes a single string, the string has N comma-separated-numbers where N is the number of instances passed"
        nargs = 1
      "--div-and-round-nearest"
        help = "divide instance lenght and width (but not profit) by the passed factor and round them to nearest (THE ALGORITHM BECOMES A GUESS, DO NOT GIVE A PRIMAL SOLUTION NOR A VALID DUAL SOLUTION)"
        arg_type = Int
        default = [1]
        nargs = 1
      "--div-and-round-up"
        help = "divide instance lenght and width (but not profit) by the passed factor and round them up (THE ALGORITHM BECOMES A PRIMAL HEURISTIC)"
        arg_type = Int
        default = [1]
        nargs = 1
      "--div-and-round-down"
        help = "divide instance lenght and width (but not profit) by the passed factor and round them down (THE ALGORITHM BECOMES A OPTIMISTIC GUESS, VALID BOUND)"
        arg_type = Int
        default = [1]
        nargs = 1
      "instfnames"
        help = "the paths to the instances to be solved"
        action = :store_arg
        nargs = '*' # can pass no instances to call "-h"
    end

    p_args = parse_args(args, s)

    for option in ["lower-bounds", "upper-bounds"]
      if !isempty(p_args[option])
        @assert isone(length(p_args[option]))
        str_list = p_args[option][1]
        if match(r"^([1-9][0-9]*|0)(,([1-9][0-9]*|0))*$", str_list) === nothing
          @error "the --$(option) does not follow the desired format: integer,integer,integer,... (no comma after last number)"
        end
        num_list = parse.(Int64, split(str_list, ','))
        if length(num_list) != length(p_args["instfnames"])
          @error "the length of the list passed to --$(option) should be the exact same size as the number of instances provided"
        end
        p_args[option] = num_list
      end
    end

    num_rounds = sum(.! isone.(map(flag -> p_args[flag][1],
      ["div-and-round-nearest", "div-and-round-up", "div-and-round-down"]
    )))
    num_rounds > 1 && @error("only one of --div-and-round-{nearest,up,down} may be passed at the same time (what the fuck you expected?)")
    is_revised_furini = !p_args["faithful2furini2016"] && !p_args["flow-model"]
    p_args["final-pricing"] && !is_revised_furini && @error(
      "the final pricing technique is implemented just for Revised Furini model as of now"
    )
    p_args["final-pricing"] && isempty(p_args["lower-bounds"]) &&
      !p_args["warm-start"] && @error(
      "the flag --final-pricing only makes sense if a lower bound is provided (either directly by --lower-bound or indirectly by --warm-start)"
    )
    p_args["final-pricing"] && p_args["relax2lp"] && @error(
      "the flags --final-pricing and --relax2lp should not be used together; it is not clear what they should do, and the best interpretation (solving the relaxed model and doing the final pricing, without solving the unrelaxed reduced model after) is not specially useful and need extra code to work that is not worth it"
    )
    p_args["det-time-limit"] != [1.0E+75] && !using_cplex && @error(
      "only CPLEX has deterministic time, unfortunately"
    )

    return p_args
end

# Parse the command line arguments, and call the solve for each instance.
function run_batch(args = ARGS)
  @show args
  p_args = parse_script_args(args)

  if !p_args["do-not-mock-first-for-compilation"] && !isempty(p_args["instfnames"])
    read_build_solve_and_print(
      p_args, 1; no_csv_out = true, no_solver_out = true
    )
  end
  for inst_idx in eachindex(p_args["instfnames"])
    total_instance_time = @elapsed read_build_solve_and_print(
      p_args, inst_idx; no_solver_out = p_args["disable-solver-output"]
    )
    @show total_instance_time
  end
end

run_batch()

