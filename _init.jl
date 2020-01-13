push!(LOAD_PATH, "./src")
using Pkg
Pkg.activate(".")
Pkg.instantiate()
using Revise, JuMP, CPLEX, GuillotineModels, GC2DInstanceReader
m = JuMP.Model(with_optimizer(CPLEX.Optimizer, CPX_PARAM_TILIM = 10))
L, W, l, w, p, d = GC2DInstanceReader.read_instance("./test/instances/CW4")

