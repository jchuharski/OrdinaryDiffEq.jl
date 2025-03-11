module DAETS_utils
using Symbolics
using SymbolicUtils
using Hungarian
using JuMP
using GLPK
using OrdinaryDiffEq: ODEFunction, DAEFunction
using SciMLBase
using LinearAlgebra: UniformScaling
using TaylorDiff: TaylorScalar, extract_derivative, derivative
using ForwardDiff

export signature_matrix, highest_value_transversal, find_offsets, system_jacobian

"""
    signature_matrix(eqs, vars, t_ind) -> Matrix{Float64}

Construct the signature matrix Σ for the system of equations eqs with respect to
the variables vars. Each entry Σ[i,j] is the (highest) derivative order of vars[j]
appearing in eqs[i], or -Inf if it does not appear.
"""
function signature_matrix(eqs, vars, t_ind)
    if vars isa Vector{<:SymbolicUtils.BasicSymbolic} && eqs isa Vector{<:SymbolicUtils.BasicSymbolic}
        n_eqs = length(eqs)
        n_vars = length(vars)
        Σ = fill(-Inf, n_eqs, n_vars)
        
        for i in 1:n_eqs
            for j in 1:n_vars
                # Check for each variable in the equation what the highest derivative order is
                order = max_derivative_order(eqs[i], vars[j], t_ind)
                Σ[i,j] = order
            end
        end
        return Σ
    else
        println("Non-symbolic inputs detected")
        # Just use a default signature matrix for now
        # TODO: Add a better way to handle this
        n = length(vars)
        Σ = zeros(n, n)
        for i in 1:n
            for j in 1:n
                if i == j
                    Σ[i,j] = 1.0  # Assume first-order derivatives on diagonal
                else
                    Σ[i,j] = 0.0  # Assume direct dependencies off-diagonal
                end
            end
        end
        return Σ
    end
end


"""
    max_derivative_order(ex, var, t_ind)

Returns the highest derivative order of `var` that appears in the symbolic expression `ex`,
using `t_ind` as the independent variable (e.g., time). If `var` does not appear, returns -Inf.
"""
function max_derivative_order(ex, var, t_ind)
    # This is also for symbolic expressions. (probably do not need this but keeping it here for now)
    # Base case: if ex is exactly var(t_ind), order is 0.
    if isequal(ex, var(t_ind))
        return 0
    end

    # If it's a number or unrelated symbol ignore
    if ex isa Number || (ex isa Symbol && ex != var)
        return -Inf
    end

    # Check if ex is a derivative expression.
    if iscall(ex) && operation(ex) isa Symbolics.Differential
        inner = arguments(ex)[1]  # Function being differentiated.
        # Recurse
        sub_order = max_derivative_order(inner, var, t_ind)
        return sub_order == -Inf ? -Inf : 1 + sub_order
    end

    # For composite expressions (e.g., sums, products), traverse arguments
    if iscall(ex)
        best = -Inf
        for arg in arguments(ex)
            # Recursively check the order of each component
            best = max(best, max_derivative_order(arg, var, t_ind))
        end
        return best
    end
    return -Inf
end

"""
    highest_value_transversal(Σ) -> (Vector{Tuple{Int, Int}}, Float64)

Finds the highest value transversal (HVT) of the signature matrix `Σ` using the Hungarian algorithm.
Returns the transversal as a vector of (row, column) indices and its value.
"""
function highest_value_transversal(Σ::Matrix{Float64})
    n = size(Σ, 1)
    
    # The Hungarian algorithm minimizes so multiply by -1 to max
    cost_matrix = -Σ
    assignment = hungarian(cost_matrix)[1]
    # Extract transversal and its value.
    transversal = [(i, assignment[i]) for i in 1:n]
    value = sum(Σ[i, assignment[i]] for i in 1:n)
    
    return transversal, value
end


"""
    find_offsets(Σ::Matrix{Float64}, T::Vector{Tuple{Int, Int}}) -> (Vector{Int}, Vector{Int})

Finds the canonical offsets `c` and `d` for the signature matrix `Σ` and the highest value transversal `T`.
Returns the vectors `c` and `d`.
"""
function find_offsets(Σ::Matrix{Float64}, T::Vector{Tuple{Int, Int}})
    n = size(Σ, 1)
    # Create JuMP model with GLPK solver
    model = Model(GLPK.Optimizer)
    # Define variables for c and d (offsets)f
    @variable(model, c[1:n] >= 0, Int)
    @variable(model, d[1:n] >= 0, Int)
    
    # Add constraints for all i, j: d[j] - c[i] >= Σ[i, j]
    for i in 1:n
        for j in 1:n
            if Σ[i, j] != -Inf
                @constraint(model, d[j] - c[i] >= Σ[i, j])
            end
        end
    end
    
    # Add constraints for equality over transversal
    for (i, j) in T
        @constraint(model, d[j] - c[i] == Σ[i, j])
    end
    
    # min sum c and d finds the "canonical" offsets
    @objective(model, Min, sum(c) + sum(d))
    optimize!(model)
    c_values = value.(c)
    d_values = value.(d)
    return c_values, d_values
end


"""
    system_jacobian(eqs, vars, t_ind, c, d, Σ, integrator=nothing) -> Matrix

Constructs the System Jacobian matrix J for the system of equations `eqs` with respect to the variables `vars`.
The offsets `c` and `d` and the signature matrix `Σ` are used to determine the structure of the Jacobian. Trying to do
this with TaylorDiff.
"""
function system_jacobian(
    eqs,
    vars,
    t_ind,
    c::Vector{Int},
    d::Vector{Int},
    Σ::Matrix{Float64},
    integrator=nothing
)
    n = size(Σ, 1)
    J = zeros(n, n)
    
    for i in 1:n
        for j in 1:n
            if d[j] - c[i] == Σ[i,j] && Σ[i,j] != -Inf
                σ_ij = Int(Σ[i,j])
                p_actual = integrator !== nothing ? integrator.p : nothing
                
                # Create a function that evaluates the i-th equation with respect to the j-th variable
                function f_ij(x)
                    x_full = copy(vars) # copy state vector
                    x_full[j] = x #Replace j with x
                    # Evaluate
                    result = zeros(n)
                    eqs.f(result, x_full, p_actual, t_ind)
                    return result[i]
                end
                # Calculate the σ_ij-th derivative
                J[i,j] = derivative(f_ij, vars[j], Val{σ_ij}())
            else
                J[i,j] = 0.0
            end
        end
    end
    println("System Jacobian: ", J)
    return J
end

end