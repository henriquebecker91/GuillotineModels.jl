using Test

#@testset "Relax/Fix Tests" begin include("run_relax_fix_tests.jl") end
#@testset "CutPattern Tests" begin include("run_cut_pattern_tests.jl") end
@testset "Model Solving Tests" begin include("run_model_solving_tests.jl") end

