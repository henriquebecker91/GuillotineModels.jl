# Check if the necessary packages are installed and define the project
# to make them available to importation by "using".
using Pkg
Pkg.activate("..")
Pkg.instantiate()

# Load the necessary packages. The ones after the `push!` are from this
# repository.
using JuMP, CPLEX, ArgParse
push!(LOAD_PATH, "../src/")
using AllSubplatesModel, GC2DInstanceReader

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
function new_empty_and_configured_model(p_args; disable_output = false)
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
      CPX_PARAM_SCRIND = disable_output ? CPLEX.CPX_OFF : CPLEX.CPX_ON,
      CPX_PARAM_THREADS = threads,
      CPX_PARAM_RANDOMSEED = seed
    ))
  else
    m = JuMP.Model(with_optimizer(CPLEX.Optimizer,
      CPX_PARAM_DETTILIM = det_time_limit,
      #CPX_PARAM_VARSEL = CPLEX.CPX_VARSEL_MAXINFEAS,
      CPX_PARAM_SCRIND = disable_output ? CPLEX.CPX_OFF : CPLEX.CPX_ON,
      CPX_PARAM_THREADS = threads,
      CPX_PARAM_RANDOMSEED = seed
    ))
  end

  return m
end

# Read the instance, build the model, solve the model, and print related stats.
function read_build_solve_and_print(
  p_args, instfname_idx; disable_output = false
)
  instfname = p_args["instfnames"][instfname_idx]
  L, W, l, w, p, d = GC2DInstanceReader.read_instance(instfname)
  before_model_build = time()
  m = new_empty_and_configured_model(p_args; disable_output = disable_output)
  #_, hvcuts, pli_lwb, np = AllSubplatesModel.build(m, d, p, l, w, L, W;
  _, hvcuts, pli2lwsb = AllSubplatesModel.build(m, d, p, l, w, L, W;
    break_hvcut_symm = p_args["break-hvcut-symmetry"],
    only_binary = p_args["only-binary-variables"]
  )
  time_to_build_model = time() - before_model_build 
  if !disable_output
    @show instfname
    @show time_to_build_model
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
  if !disable_output
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
    if p_args["break-hvcut-symmetry"]
      ps = m[:pieces_sold]
      cm = m[:cuts_made]
      ps_nz = [(i, value(ps[i])) for i = 1:length(ps) if value(ps[i]) > 0.001]
      cm_nz = [(i, value(cm[i])) for i = 1:length(cm) if value(cm[i]) > 0.001]
      @show ps_nz
      @show cm_nz
      foreach(ps_nz) do e
        pii, _ = e
        println((l[pii], w[pii]) => (pii, d[pii], p[pii]))
      end
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
      ps = m[:picuts]
      cm = m[:cuts_made]
      ps_nz = [(i, value(ps[i])) for i = 1:length(ps) if value(ps[i]) > 0.001]
      cm_nz = [(i, value(cm[i])) for i = 1:length(cm) if value(cm[i]) > 0.001]
      @show ps_nz
      @show cm_nz
      foreach(ps_nz) do e
        i, _ = e
        pli, pii = np[i]
        println((pli_lwb[pli][1], pli_lwb[pli][2]) => (l[pii], w[pii]))
      end
      foreach(cm_nz) do e
        i, _ = e
        parent, fchild, schild = hvcuts[i]
        println((pli_lwb[parent][1], pli_lwb[parent][2]) => ((pli_lwb[fchild][1], pli_lwb[fchild][2]), (pli_lwb[schild][1], pli_lwb[schild][2])))
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
      "--do-not-mock-first-for-compilation"
        help = "if passed then the first instance will not be solved twice"
        nargs = 0
      "--break-hvcut-symmetry"
        help = "will double the number of variables, and reduce symmetry"
        nargs = 0
      "--only-binary-variables"
        help = "CAUTION DO NOT GUARANTEE OPTIMALITY FOR NOW, IS HEURISTIC"
        nargs = 0
      "instfnames"
        help = "the paths to the instances to be solved"
        action = :store_arg
        nargs = '*' # can pass no instances to call "-h"
    end

    return parse_args(args, s)
end

# Parse the command line arguments, and call the solve for each instance.
function run_batch(args = ARGS)
  @show args
  p_args = parse_script_args(args)

  if !p_args["do-not-mock-first-for-compilation"] && !isempty(p_args["instfnames"])
    read_build_solve_and_print(p_args, 1, disable_output = true)
  end
  for inst_idx in eachindex(p_args["instfnames"])
    read_build_solve_and_print(p_args, inst_idx)
  end
end

run_batch()

