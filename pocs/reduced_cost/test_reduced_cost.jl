using JuMP, GLPK, Test, MathOptInterface
const MOI = MathOptInterface

function reduced_cost(x::VariableRef)::Float64
    !has_duals(owner_model(x)) && error(
        "The reduced cost is not available because no dual result is" *
        " available."
    )
    is_fixed(x) && return shadow_price(FixRef(x))
    rc = 0.0
    has_upper_bound(x) && (rc += shadow_price(UpperBoundRef(x)))
    has_lower_bound(x) && (rc += shadow_price(LowerBoundRef(x)))
    return rc
end

function Base.empty!(model::Model)::Model
    # The method changes the Model object to, basically, the state it was when
    # created (if the optimizer was already pre-configured). The two exceptions
    # are: optimize_hook and operator_counter. One is left alone because it is
    # related to optimizer attributes and the other is just a counter for a
    # single-time warning message (so keeping it helps to discover
    # inneficiencies).
    MOI.empty!(model.moi_backend)
    empty!(model.shapes)
    empty!(model.bridge_types)
    model.nlp_data = nothing
    empty!(model.obj_dict)
    empty!(model.ext)
    return model
end

function pretty_print_model(model)
    iob = IOBuffer()
    write(iob, model)
    model_str = String(take!(iob))
    println(model_str)
end

function test_reduced_costs()
    @testset "test_reduced_costs" begin
    model = JuMP.Model(GLPK.Optimizer)
    let
        @variable(model, a >= 0)
        @objective(model, Min, a)
				#pretty_print_model(model)
				@test_throws ErrorException reduced_cost(a)
				optimize!(model)
				@test reduced_cost(a) ≈ -1.0 atol=1e-6 rtol=1e-6
    end
    Base.empty!(model)
    let
        @variable(model, a >= 0)
        @objective(model, Min, 0.5 * a)
				#pretty_print_model(model)
				@test_throws ErrorException reduced_cost(a)
				optimize!(model)
				@test reduced_cost(a) ≈ -0.5 atol=1e-6 rtol=1e-6
    end
    Base.empty!(model)
    let
        @variable(model, a >= 0)
        @objective(model, Max, -1 * a)
				#pretty_print_model(model)
				@test_throws ErrorException reduced_cost(a)
				optimize!(model)
				@test reduced_cost(a) ≈ 1 atol=1e-6 rtol=1e-6
    end
    Base.empty!(model)
    let
        @variable(model, 10 >= a >= 0)
        @objective(model, Max, -1 * a)
				#pretty_print_model(model)
				@test_throws ErrorException reduced_cost(a)
				optimize!(model)
				@test reduced_cost(a) ≈ 1 atol=1e-6 rtol=1e-6
    end
    Base.empty!(model)
    let
        @variable(model, 0 <= a <= 10)
        @objective(model, Max, a)
				#pretty_print_model(model)
				@test_throws ErrorException reduced_cost(a)
				optimize!(model)
				@show value(a)
				@show reduced_cost(a)
				@show shadow_price(LowerBoundRef(a))
				@show shadow_price(UpperBoundRef(a))
				#@test reduced_cost(a) ≈ 0.0 atol=1e-6 rtol=1e-6
    end
		exit()
    Base.empty!(model)
    let
        @variable(model, a >= 0)
        @variable(model, b >= 0)
				@constraint(model, a + b <= 5)
        @objective(model, Max, 2a + b)
				#pretty_print_model(model)
				@test_throws ErrorException reduced_cost(a)
				optimize!(model)
				@test reduced_cost(a) ≈ 0.0 atol=1e-6 rtol=1e-6
				@test reduced_cost(b) ≈ 1.0 atol=1e-6 rtol=1e-6
    end
    Base.empty!(model)
    let
        @variable(model, a >= 0)
        @variable(model, b >= 0)
				@constraint(model, a + b >= -5)
        @objective(model, Min, -2a + -b)
				#pretty_print_model(model)
				@test_throws ErrorException reduced_cost(a)
				optimize!(model)
				@test reduced_cost(a) ≈ 0.0 atol=1e-6 rtol=1e-6
				@test reduced_cost(b) ≈ -1.0 atol=1e-6 rtol=1e-6
    end
    Base.empty!(model)
    end
end

test_reduced_costs()

