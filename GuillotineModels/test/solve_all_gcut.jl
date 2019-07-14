using JuMP, CPLEX
push!(LOAD_PATH, "../src")
using AllSubplatesModel, GC2DInstanceReader

#=
fname = "./instances/okp2"
@show fname
m = JuMP.Model(with_optimizer(CPLEX.Optimizer,
  #CPX_PARAM_TILIM = 300,
  #CPX_PARAM_PROBE = -1,
  CPX_PARAM_THREADS = 1#,
#  CPX_PARAM_RANDOMSEED = 0
))
L, W, l, w, p, d = GC2DInstanceReader.read_instance(fname)
_, ilwb, nnn, np = AllSubplatesModel.build(m, d, p, l, w, L, W)
@show m
optimize!(m)
@show JuMP.objective_value(m)
=#
#=
fname = "./instances/CHL7"
@show fname
m = JuMP.Model(with_optimizer(CPLEX.Optimizer,
  #CPX_PARAM_TILIM = 300,
  #CPX_PARAM_PROBE = -1,
  CPX_PARAM_THREADS = 1,
  CPX_PARAM_RANDOMSEED = 0
))
L, W, l, w, p, d = GC2DInstanceReader.read_instance(fname)
p = l .* w # NOTE: this is done because CHL7 does not include profits
_, ilwb, nnn, np = AllSubplatesModel.build(m, d, p, l, w, L, W)
@show m
optimize!(m)
@show JuMP.objective_value(m)
=#
for i = 1:13
  fname = "./instances/gcut$(i).txt"
  @show fname
  m = JuMP.Model(with_optimizer(CPLEX.Optimizer,
    CPX_PARAM_THREADS = 1,
    CPX_PARAM_RANDOMSEED = 0
  ))
  L, W, l, w, p, d = GC2DInstanceReader.read_instance(fname)
  AllSubplatesModel.build(m, d, p, l, w, L, W)
  @show m
  optimize!(m)
  @show JuMP.objective_value(m)
end

