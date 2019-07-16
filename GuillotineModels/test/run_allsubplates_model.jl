using JuMP, CPLEX, ArgParse
push!(LOAD_PATH, "../src")
using AllSubplatesModel, GC2DInstanceReader

function new_empty_and_configured_model(p_args; disable_output = false)
  @assert isone(length(p_args["threads"]))
  threads = first(p_args["threads"])
  @assert isone(length(p_args["seed"]))
  seed = first(p_args["seed"])

  if p_args["disable-cplex-tricks"]
    m = JuMP.Model(with_optimizer(CPLEX.Optimizer,
      CPX_PARAM_PROBE = -1,
      CPX_PARAM_HEURFREQ = -1,
      CPX_PARAM_REPEATPRESOLVE = -1,
      CPX_PARAM_DIVETYPE = 1,
      CPX_PARAM_SCRIND = disable_output ? CPLEX.CPX_OFF : CPLEX.CPX_ON,
      CPX_PARAM_THREADS = threads,
      CPX_PARAM_RANDOMSEED = seed
    ))
  else
    m = JuMP.Model(with_optimizer(CPLEX.Optimizer,
      CPX_PARAM_SCRIND = disable_output ? CPLEX.CPX_OFF : CPLEX.CPX_ON,
      CPX_PARAM_THREADS = threads,
      CPX_PARAM_RANDOMSEED = seed
    ))
  end

  return m
end

function read_build_solve_and_print(
  p_args, instfname_idx; disable_output = false
)
  instfname = p_args["instfnames"][instfname_idx]
  L, W, l, w, p, d = GC2DInstanceReader.read_instance(instfname)
  before_model_build = time()
  m = new_empty_and_configured_model(p_args; disable_output = disable_output)
  _, ilwb, nnn, np = AllSubplatesModel.build(m, d, p, l, w, L, W)
  time_to_build_model = time() - before_model_build 
  !disable_output && @show instfname
  !disable_output && @show time_to_build_model
  time_to_solve_model = @elapsed optimize!(m)
  !disable_output && @show time_to_solve_model
  !disable_output && @show objective_value(m)

  return nothing
end

function parse_script_args(args = ARGS)
    s = ArgParseSettings()
    @add_arg_table s begin
      "--threads"
        help = "number of threads used by CPLEX (not model building)"
        arg_type = Int
        default = [1]
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
      "instfnames"
        help = "the paths to the instances to be solved"
        action = :store_arg
        nargs = '+'
        required = true
    end

    return parse_args(args, s)
end

function run_batch(args = ARGS)
  @show args
  p_args = parse_script_args(args)

  if !p_args["do-not-mock-first-for-compilation"] 
    read_build_solve_and_print(p_args, 1, disable_output = true)
  end
  for inst_idx in eachindex(p_args["instfnames"])
    read_build_solve_and_print(p_args, inst_idx)
  end
end

run_batch()

