# Check if the necessary packages are installed and define the project
# to make them available to importation by "using".
using Pkg
Pkg.activate("..")
Pkg.instantiate()

# Load the necessary packages. The ones after the `push!` are from this
# repository.
using JuMP, CPLEX, ArgParse, MathOptInterface
push!(LOAD_PATH, "../src/")
using AllSubplatesModel, FlowModel, GC2DInstanceReader

# JuMP has no method for getting all constraints, you need to get the
# types of constraints used in the model and then query the number for
# specific type.
function num_all_constraints(m) :: Int64
  sum = 0 :: Int64
  for (ftype, stype) in list_of_constraint_types(m)
    sum += num_constraints(m, ftype, stype)
  end
  return sum
end

# Create a new CPLEX model and set its configurations.
function new_empty_and_configured_model(p_args; no_cplex_out = no_cplex_out)
  @assert isone(length(p_args["threads"]))
  threads = first(p_args["threads"])
  @assert isone(length(p_args["seed"]))
  seed = first(p_args["seed"])
  @assert isone(length(p_args["det-time-limit"]))
  det_time_limit = first(p_args["det-time-limit"])

  if p_args["disable-cplex-tricks"]
    m = JuMP.Model(with_optimizer(CPLEX.Optimizer,
      CPX_PARAM_PROBE = -1,
      CPX_PARAM_HEURFREQ = -1,
      CPX_PARAM_REPEATPRESOLVE = -1,
      CPX_PARAM_DIVETYPE = 1,
      CPX_PARAM_DETTILIM = det_time_limit,
      #CPX_PARAM_VARSEL = CPLEX.CPX_VARSEL_MAXINFEAS,
      CPX_PARAM_SCRIND = no_cplex_out ? CPLEX.CPX_OFF : CPLEX.CPX_ON,
      CPX_PARAM_THREADS = threads,
      CPX_PARAM_RANDOMSEED = seed
    ))
  else
    m = JuMP.Model(with_optimizer(CPLEX.Optimizer,
      CPX_PARAM_DETTILIM = det_time_limit,
      #CPX_PARAM_VARSEL = CPLEX.CPX_VARSEL_MAXINFEAS,
      CPX_PARAM_SCRIND = no_cplex_out ? CPLEX.CPX_OFF : CPLEX.CPX_ON,
      CPX_PARAM_THREADS = threads,
      CPX_PARAM_RANDOMSEED = seed
    ))
  end

  return m
end

# Read the instance, build the model, solve the model, and print related stats.
function read_build_solve_and_print(
  p_args, instfname_idx; no_csv_out = false, no_cplex_out = false
)
  instfname = p_args["instfnames"][instfname_idx]
  L, W, l, w, p, d = GC2DInstanceReader.read_instance(instfname)
  before_model_build = time()
  m = new_empty_and_configured_model(p_args; no_cplex_out = no_cplex_out)
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
        lb = get(p_args["lower-bounds"], instfname_idx, 0),
        ub = get(p_args["upper-bounds"], instfname_idx, sum(d .* p))
      )
    end
  end
  time_to_build_model = time() - before_model_build 
  if !no_csv_out
    @show instfname
    seed = first(p_args["seed"])
    @show seed
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
  end
  flush(stdout) # guarantee all messages will be flushed before calling cplex
  time_to_solve_model = @elapsed optimize!(m)
  sleep(0.01) # shamefully needed to avoid out of order messages from cplex
  if !no_csv_out
    @show time_to_solve_model
    if primal_status(m) == MOI.FEASIBLE_POINT
      obj_value = objective_value(m)
    else
      obj_value = 0
    end
    @show obj_value
    obj_bound = objective_bound(m)
    @show obj_bound
    stop_reason = termination_status(m)
    @show stop_reason
    stop_code = Int64(stop_reason)
    @show stop_code
    if p_args["flow-model"]
      ps = value.(m[:edge][1:length(d)])
      ps_nz = [iv for iv in enumerate(ps) if iv[2] > 0.001]
      @show ps_nz
    else
      if p_args["break-hvcut-symmetry"]
        ps = m[:pieces_sold]
        cm = m[:cuts_made]
        ps_nz = [(i, value(ps[i])) for i = 1:length(ps) if value(ps[i]) > 0.001]
        cm_nz = [(i, value(cm[i])) for i = 1:length(cm) if value(cm[i]) > 0.001]
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
      else
        num_plates = length(pli_lwb)
        @show num_plates
        ps = m[:picuts]
        num_picuts = length(ps)
        @show num_picuts
        cm = m[:cuts_made]
        num_cuts_made = length(cm)
        @show num_cuts_made
        ps_nz = [(i, value(ps[i])) for i = 1:length(ps) if value(ps[i]) > 0.001]
        cm_nz = [(i, value(cm[i])) for i = 1:length(cm) if value(cm[i]) > 0.001]
        @show ps_nz
        @show cm_nz
        println("(plate length, plate width) => (number of times this extraction happened, piece length, piece width)")
        foreach(ps_nz) do e
          i, v = e
          pli, pii = np[i]
          println((pli_lwb[pli][1], pli_lwb[pli][2]) => (v, l[pii], w[pii]))
        end
        println("(parent plate length, parent plate width) => ((first child plate length, first child plate width), (second child plate length, second child plate width))")
        foreach(cm_nz) do e
          i, _ = e
          parent, fchild, schild = hvcuts[i]
          println((pli_lwb[parent][1], pli_lwb[parent][2]) => ((pli_lwb[fchild][1], pli_lwb[fchild][2]), (pli_lwb[schild][1], pli_lwb[schild][2])))
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
        help = "number of threads used by CPLEX (not model building)"
        arg_type = Int
        default = [1]
        nargs = 1
      "--det-time-limit"
        help = "deterministic time limit for CPLEX (not model building)"
        arg_type = Float64
        default = [1.0E+75]
        nargs = 1
      "--seed"
        help = "seed used by CPLEX (model building is deterministic)"
        arg_type = Int
        default = [1]
        nargs = 1
      "--disable-cplex-tricks"
        help = "disable heuristics, probe, repeat presolve and dive"
        nargs = 0
      "--disable-cplex-output"
        help = "disable the CPLEX output for all instances"
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

    return p_args
end

# Parse the command line arguments, and call the solve for each instance.
function run_batch(args = ARGS)
  @show args
  p_args = parse_script_args(args)

  if !p_args["do-not-mock-first-for-compilation"] && !isempty(p_args["instfnames"])
    read_build_solve_and_print(
      p_args, 1; no_csv_out = true, no_cplex_out = true
    )
  end
  for inst_idx in eachindex(p_args["instfnames"])
    read_build_solve_and_print(
      p_args, inst_idx; no_cplex_out = p_args["disable-cplex-output"]
    )
  end
end

run_batch()

