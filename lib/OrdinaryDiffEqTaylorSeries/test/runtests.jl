using OrdinaryDiffEqTaylorSeries, ODEProblemLibrary, DiffEqDevTools
using Test
using Symbolics
using Plots
using LinearAlgebra  # Add this import for the norm function

# @testset "Taylor2 Convergence Tests" begin
#     # Test convergence
#     dts = 2. .^ (-8:-4)
#     testTol = 0.2
#     sim = test_convergence(dts, prob_ode_linear, ExplicitTaylor2())
#     @test sim.𝒪est[:final]≈2 atol=testTol
#     sim = test_convergence(dts, prob_ode_2Dlinear, ExplicitTaylor2())
#     @test sim.𝒪est[:final]≈2 atol=testTol
# end

# @testset "TaylorN Convergence Tests" begin
#     # Test convergence
#     dts = 2. .^ (-8:-4)
#     testTol = 0.2
#     for N in 3:4
#         alg = ExplicitTaylor(order=Val(N))
#         sim = test_convergence(dts, prob_ode_linear, alg)
#         @test sim.𝒪est[:final]≈N atol=testTol
#         sim = test_convergence(dts, prob_ode_2Dlinear, alg)
#         @test sim.𝒪est[:final]≈N atol=testTol
#     end
# end

# println("DONE with ODE tests")
# # include(joinpath(@__DIR__, "../src/DAETS_utils.jl"))
# # include(joinpath(@__DIR__, "../src/TaylorSeries_caches.jl"))
# println("Starting tests on DAETS")

# # --- Test Cases --- #
# @testset "Signature Matrix & Derivative Order Tests" begin
#     @syms t x(t) y(t)
#     # # Test max_derivative_order
#     # The variable itself.
#     order1 = max_derivative_order(x(t), x, t)
#     @test order1 == 0

#     # A multiplication expression: 7*x(t)
#     order2 = max_derivative_order(7*x(t), x, t)
#     @test order2 == 0

#     # First derivative: x'(t)
#     dx = Differential(t)(x(t))
#     # println(typeof(dx))
#     order3 = max_derivative_order(dx, x, t)
#     @test order3 == 1

#     # Second derivative: x''(t)
#     d2x = Differential(t)(dx)
#     order4 = max_derivative_order(d2x, x, t)
#     @test order4 == 2

#     # Expression that does not contain x(t): y(t)
#     order5 = max_derivative_order(y(t), x, t)
#     @test order5 == -Inf

#     # Test signature_matrix:
#     # Equation 1: f₁ = x'(t) + sin(x(t)) = 0
#     # Equation 2: f₂ = y''(t) - x(t) = 0
#     eq1 = Differential(t)(x(t)) + sin(x(t))
#     eq2 = Differential(t)(Differential(t)(y(t))) - x(t)
#     eqs = [eq1, eq2]
#     # pass functions x and y (which take input t)
#     vars = [x, y]

#     Σ = signature_matrix(eqs, vars, t)
#     # First equation:
#     # Column 1: x(t) appears as x'(t) ⟹ order 1.
#     # Column 2: no y(t) ⟹ -Inf.
#     @test Σ[1, 1] == 1
#     @test Σ[1, 2] == -Inf

#     # Second equation:
#     # Column 1: x(t) appears as x(t) ⟹ order 0.
#     # Column 2: y(t) appears as y''(t) ⟹ order 2.
#     @test Σ[2, 1] == 0
#     @test Σ[2, 2] == 2
# end

# println("DONE Signature Matrix Tests")

# @testset "Highest Value Transversal Tests" begin
#     @syms t x(t) y(t)

#     # Same Equation
#     eq1 = Differential(t)(x(t)) + sin(x(t))
#     eq2 = Differential(t)(Differential(t)(y(t))) - x(t)
#     eqs = [eq1, eq2]
#     vars = [x, y]
#     Σ = signature_matrix(eqs, vars, t)

#     # Expected signature matrix:
#     # [ 1  -Inf
#     #   0     2 ]
#     @test Σ == [1 -Inf; 0 2]
#     # Test HVT
#     transversal, value = highest_value_transversal(Σ)

#     # Expected transversal: [(1, 1), (2, 2)] (x'(t) in eq1 and y''(t) in eq2)
#     # Expected value: 1 + 2 = 3 TODO: Should probably add a check in case the value is infinite. In this case we say that the system is "ill-posed".
#     @test transversal == [(1, 1), (2, 2)]
#     @test value == 3
# end
# @testset "Highest Value Transversal Tests for 3x3 System" begin
#     @syms t x(t) y(t) z(t)

#     # Equation 1: f₁ = x''(t) + y(t) = 0
#     # Equation 2: f₂ = y'(t) + z(t) = 0
#     # Equation 3: f₃ = z''(t) + x(t) = 0
#     eq1 = Differential(t)(Differential(t)(x(t))) + y(t)
#     eq2 = Differential(t)(y(t)) + z(t)
#     eq3 = Differential(t)(Differential(t)(z(t))) + x(t)
#     eqs = [eq1, eq2, eq3]
#     vars = [x, y, z]
#     Σ = signature_matrix(eqs, vars, t)

#     # Expected signature matrix:
#     # [ 2   0  -Inf
#     #  -Inf  1     0
#     #   0  -Inf    2 ]
#     @test Σ == [2 0 -Inf; -Inf 1 0; 0 -Inf 2]
#     # Test HVT
#     transversal, value = highest_value_transversal(Σ)

#     # Expected transversal: [(1, 1), (2, 2), (3, 3)] (x''(t) in eq1, y'(t) in eq2, z''(t) in eq3)
#     # Expected value: 2 + 1 + 2 = 5
#     @test transversal == [(1, 1), (2, 2), (3, 3)]
#     @test value == 5
# end
# println("DONE Highest Value Transversal Tests")

# @testset "Find Offsets Tests for 2x2 System" begin
#     @syms t x(t) y(t)

#     # Same Equation as 2x2 System
#     eq1 = Differential(t)(x(t)) + sin(x(t))
#     eq2 = Differential(t)(Differential(t)(y(t))) - x(t)
#     eqs = [eq1, eq2]
#     vars = [x, y]
#     Σ = signature_matrix(eqs, vars, t)

#     # Find the highest value transversal
#     transversal, value = highest_value_transversal(Σ)

#     # Find the offsets
#     c, d = find_offsets(Σ, transversal)

#     # Expected offsets (canonical):
#     # c = [0, 0]
#     # d = [1, 2]
#     @test c == [0, 0]
#     @test d == [1, 2]
# end
# @testset "Find Offsets Tests for 3x3 System" begin
#     @syms t x(t) y(t) z(t)

#     # Same Equation as 3x3 System
#     eq1 = Differential(t)(Differential(t)(x(t))) + y(t)
#     eq2 = Differential(t)(y(t)) + z(t)
#     eq3 = Differential(t)(Differential(t)(z(t))) + x(t)
#     eqs = [eq1, eq2, eq3]
#     vars = [x, y, z]
#     Σ = signature_matrix(eqs, vars, t)
#     transversal, value = highest_value_transversal(Σ)

#     # Test Offsets
#     c, d = find_offsets(Σ, transversal)

#     # Expected offsets (canonical):
#     # c = [0, 0, 0]
#     # d = [2, 1, 2]
#     @test c == [0, 0, 0]
#     @test d == [2, 1, 2]
# end

# println("DONE Find Offsets Tests")

# @testset "System Jacobian Tests for 2x2 System" begin
#     @syms t x(t) y(t)

#     # Same 2x2
#     eq1 = Differential(t)(x(t)) + sin(x(t))
#     eq2 = Differential(t)(Differential(t)(y(t))) - x(t)
#     eqs = [eq1, eq2]
#     vars = [x, y]
#     Σ = signature_matrix(eqs, vars, t)
#     transversal, value = highest_value_transversal(Σ)
#     c, d = find_offsets(Σ, transversal)

#     # Convert c and d to Vector{Int}
#     c = Int.(c)
#     d = Int.(d)

#     # Test Jacobian
#     J = system_jacobian(eqs, vars, t, c, d, Σ)

#     # Expected Jacobian:
#     # [ 1   0
#     #   0   1 ]
#     @test isequal(J[1, 1], 1)
#     @test isequal(J[1, 2], 0)
#     @test isequal(J[2, 1], 0)
#     @test isequal(J[2, 2], 1)
# end
# @testset "System Jacobian Tests for 3x3 System" begin
#     @syms t x(t) y(t) z(t)

#     # Same 3x3
#     eq1 = Differential(t)(Differential(t)(x(t))) + y(t)
#     eq2 = Differential(t)(y(t)) + z(t)
#     eq3 = Differential(t)(Differential(t)(z(t))) + x(t)
#     eqs = [eq1, eq2, eq3]
#     vars = [x, y, z]
#     Σ = signature_matrix(eqs, vars, t)
#     transversal, value = highest_value_transversal(Σ)
#     c, d = find_offsets(Σ, transversal)

#     # Convert c and d to Vector{Int}
#     c = Int.(c)
#     d = Int.(d)

#     # Test Jacobian
#     J = system_jacobian(eqs, vars, t, c, d, Σ)

#     # Expected Jacobian:
#     # [1   0   0
#     #  0   1   1
#     #  0   0   1]
#     @test isequal(J[1, 1], 1)
#     @test isequal(J[1, 2], 0)
#     @test isequal(J[1, 3], 0)
#     @test isequal(J[2, 1], 0)
#     @test isequal(J[2, 2], 1)
#     @test isequal(J[2, 3], 0)
#     @test isequal(J[3, 1], 0)
#     @test isequal(J[3, 2], 0)
#     @test isequal(J[3, 3], 1)
# end
# println("DONE System Jacobian Tests")

# @testset "System Jacobian Tests for Simple Pendulum" begin
#     @syms t x(t) y(t) λ(t) G L

#     # Pendulum Equations
#     f = Differential(t)(Differential(t)(x(t))) + x(t) * λ(t)
#     g = Differential(t)(Differential(t)(y(t))) + y(t) * λ(t) - G
#     h = x(t)^2 + y(t)^2 - L^2
#     eqs = [f, g, h]
#     vars = [x, y, λ]
#     Σ = signature_matrix(eqs, vars, t)
#     transversal, value = highest_value_transversal(Σ)
#     c, d = find_offsets(Σ, transversal)
#     c = Int.(c)
#     d = Int.(d)
#     J = system_jacobian(eqs, vars, t, c, d, Σ)

#     # Expected Jacobian:
#     # [1   0   x(t)
#     #  0   1   y(t)
#     #  2   2   0]
#     @test isequal(J[1, 1], 1)
#     @test isequal(J[1, 2], 0)
#     @test isequal(J[1, 3], x(t))
#     @test isequal(J[2, 1], 0)
#     @test isequal(J[2, 2], 1)
#     @test isequal(J[2, 3], y(t))
#     @test isequal(J[3, 1], 2x(t))
#     @test isequal(J[3, 2], 2y(t))
#     @test isequal(J[3, 3], 0)
# end

using OrdinaryDiffEq
using Test
function exponential_decay!(du, u, p, t)
    du[1] = -u[1]  # du/dt = -u
end

# Initial condition
u0 = [1.0]  # u(0) = 1
tspan = (0.0, 1.0)  # Time span from 0 to 1
prob = ODEProblem(exponential_decay!, u0, tspan)
dt = 0.01  # Fixed timestep
sol = solve(prob, DAETS(), dt=dt)
exact_solution(t) = exp(-t)
for t in sol.t
    u_num = sol(t)[1]
    u_exact = exact_solution(t)
    @test abs(u_num - u_exact) < 5  # Verify accuracy (up to first order)
    println("total_error: ", abs(u_num - u_exact))
end


####### Expected outputs ######
### Signature Matrix: [1]
### Jacobian Matrix: [1]
# Create a plot comparing solutions
t_dense = range(tspan[1], tspan[2], length=100)
exact_values = exact_solution.(t_dense)
p = plot(t_dense, exact_values, label="Exact Solution", linewidth=2, 
         title="Exponential Decay: Numerical vs Exact Solution",
         xlabel="Time", ylabel="u(t)", legend=:topright)
         
plot!(p, sol.t, [u[1] for u in sol.u], label="DAETS Solution", 
      marker=:circle, markersize=3, linestyle=:dash)
savefig(p, joinpath(@__DIR__, "exponential_decay_comparison2.png"))

println("ODE test passed!")

# Damped Harmonic Oscillator

####### Expected outputs ######
### Signature Matrix: [1, 0; 0, 1]
### Jacobian Matrix: [1, 0; 0, 1]

@testset "Damped Harmonic Oscillator" begin
    using LinearAlgebra
    function harmonic_oscillator!(du, u, p, t)
        ω₀, ζ = p  # Natural frequency and damping ratio
        x, v = u   # Position and velocity
        
        du[1] = v                  # x' = v
        du[2] = -ω₀^2 * x - 2*ζ*ω₀*v  # v' = -ω₀²x - 2ζω₀v
    end
    p = [2π, 0.1]
    u0 = [1.0, 0.0]
    tspan = (0.0, 5.0)
    
    prob = ODEProblem(harmonic_oscillator!, u0, tspan, p)
    
    # Solve with DAETS
    dt = 0.01
    sol = solve(prob, DAETS(), dt=dt)
    
    # Analytical solution for underdamped case (ζ < 1)
    function exact_solution(t, ω₀, ζ)
        # Damped frequency
        ωd = ω₀ * sqrt(1 - ζ^2)
        
        # Solution for x(t) with x(0) = 1, v(0) = 0
        x = exp(-ζ*ω₀*t) * (cos(ωd*t) + (ζ*ω₀/ωd)*sin(ωd*t))
        
        # Solution for v(t)
        v = -ω₀ * exp(-ζ*ω₀*t) * (sin(ωd*t) + ζ*cos(ωd*t))
        
        return [x, v]
    end

    test_times = range(tspan[1], tspan[2], length=10)
    max_error = 0.0
    
    for t in test_times
        u_num = sol(t)
        u_exact = exact_solution(t, p[1], p[2])
        
        error = norm(u_num - u_exact)
        max_error = max(max_error, error)
        
        println("t = $t, Error: $error")
    end
    @test max_error < 10
    
    # Create a plot comparing solutions
    t_dense = range(tspan[1], tspan[2], length=200)
    exact_positions = [exact_solution(t, p[1], p[2])[1] for t in t_dense]
    
    p = plot(t_dense, exact_positions, label="Exact Solution", linewidth=2, 
             title="Damped Harmonic Oscillator: Position vs Time",
             xlabel="Time (s)", ylabel="Position", legend=:topright)
             
    plot!(p, sol.t, [u[1] for u in sol.u], label="DAETS Solution", 
          marker=:circle, markersize=2, linestyle=:dash)
    
    savefig(p, joinpath(@__DIR__, "harmonic_oscillator2.png"))
    
    println("Harmonic oscillator test passed!")
end

