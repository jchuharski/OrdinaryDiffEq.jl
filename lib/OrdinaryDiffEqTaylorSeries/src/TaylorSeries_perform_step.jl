using TaylorDiff: TaylorDiff, extract_derivative, extract_derivative!

@inline make_taylor(all::Vararg{X, P}) where {P, X <: AbstractArray} = TaylorArray(Base.first(all), Base.tail(all))
@inline make_taylor(all::Vararg{X, P}) where {P, X} = TaylorScalar(all)

function initialize!(integrator, cache::ExplicitTaylor2ConstantCache)
    integrator.kshortsize = 3
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
end

@muladd function perform_step!(integrator, cache::ExplicitTaylor2ConstantCache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    k1 = f(uprev, p, t)
    u1 = make_taylor(uprev, k1)
    t1 = TaylorScalar{1}(t, one(t))
    k2 = f(u1, p, t1).partials[1]
    u = @.. uprev + dt * k1 + dt^2 / 2 * k2
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 3)
    integrator.k[1] = k1
    integrator.k[2] = k2
    integrator.u = u
end

function initialize!(integrator, cache::ExplicitTaylor2Cache)
    integrator.kshortsize = 3
    resize!(integrator.k, integrator.kshortsize)
    # Setup k pointers
    integrator.k[1] = cache.k1
    integrator.k[2] = cache.k2
    integrator.k[3] = cache.k3
    return nothing
end

@muladd function perform_step!(integrator, cache::ExplicitTaylor2Cache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack k1, k2, k3, utilde, tmp = cache

    # The following code is written to be fully non-allocating
    f(k1, uprev, p, t)
    u1 = make_taylor(uprev, k1)
    t1 = TaylorScalar{1}(t, one(t))
    out1 = make_taylor(k1, k2)
    f(out1, u1, p, t1)
    @.. u = uprev + dt * k1 + dt^2 / 2 * k2
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 3)
    return nothing
end

function initialize!(integrator, cache::ExplicitTaylorConstantCache{P}) where P
    integrator.kshortsize = P
    integrator.k = typeof(integrator.k)(undef, P)
end

@muladd function perform_step!(integrator, cache::ExplicitTaylorConstantCache{P}, repeat_step = false) where P
    @unpack t, dt, uprev, u, f, p = integrator
    us = typeof(u)[]
    integrator.k[1] = f(uprev, p, t)
    push!(us, integrator.k[1])
    u = @.. uprev + dt * us[1]
    dti = dt
    for i in 1:P-1
        ui = make_taylor(uprev, us...)
        ti = TaylorScalar{i}(t, one(t))
        integrator.k[i + 1] = f(ui, p, ti).partials[i]
        push!(us, integrator.k[i + 1] / (i + 1))
        dti *= dt
        u += dti * us[i + 1]
    end
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, P)
    integrator.u = u
end

function initialize!(integrator, cache::ExplicitTaylorCache{P}) where P
    integrator.kshortsize = P
    resize!(integrator.k, P)
    # Setup k pointers
    for (i, k) in enumerate(cache.ks)
        integrator.k[i] = k
    end
    return nothing
end

@muladd function perform_step!(integrator, cache::ExplicitTaylorCache{P}, repeat_step = false) where P
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack ks, us = cache

    # The following code is written to be fully non-allocating
    f(ks[1], uprev, p, t)
    @.. us[1] .= ks[1]
    @.. u = uprev + dt * us[1]
    dti = dt
    for i in 1:P-1
        ui = make_taylor(uprev, us[1:i]...)
        ti = TaylorScalar{i}(t, one(t))
        outi = make_taylor(ks[1:i+1]...)
        f(outi, ui, p, ti)
        us[i + 1] .= ks[i + 1] / (i + 1)
        dti *= dt
        @.. u += dti * us[i + 1]
    end
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 3)
    return nothing
end

# DAETS
function initialize!(integrator, cache::DAETSCache)
    resize!(cache.xTS, length(integrator.u))
    fill!(cache.xTS, 0.0)
    cache.xtrial = zero(integrator.u)  # Projected solution
    cache.htrial = integrator.dt  # Step size
    cache.e = 0.0  # Error estimate
    cache.tmp = zero(integrator.u)  # Temp storage for intermediate calculations
    cache.atmp = zero(integrator.u)  # Temp storage for error compute
    f = integrator.f 
    
    # Use the state vector as variables
    vars = integrator.u
    t = integrator.t 
    # TODO: Fix this or sotre it or something we use this multiple times...
    order = 5
    if hasproperty(integrator.alg, :order)
        order = integrator.alg.order
    end
    # Preprocessing
    cache.Σ = num_signature_matrix(integrator, order)
    # cache.Σ = signature_matrix(f, vars, t)
    println("Signature matrix: ", cache.Σ)
    T, _ = highest_value_transversal(cache.Σ)
    cache.c, cache.d = find_offsets(cache.Σ, T)
    cache.J = system_jacobian(f, vars, t, cache.c, cache.d, cache.Σ, integrator)
    cache.coeffs = [zero(integrator.u) for _ in 1:order]
    return nothing
end


# Below are more helper functions for DAETS. Couldn't get these to work with the tests in a different file so I put them here. Probably should be moved to DAETS_utils.jl somehow.
function compute_taylor_coefficients!(integrator, cache)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack Σ, c, d, J, xTS, xtrial, htrial, e, tmp = cache

    order = 5  # TODO: Fix this at some point
    if hasproperty(integrator.alg, :order)
        order = integrator.alg.order
    end

    n = length(u)
    # Initialize all coefficients
    coeffs = [zeros(n) for _ in 1:order]
    
    # First coefficient is the function
    f(tmp, uprev, p, t)
    coeffs[1] .= tmp
    
    # Higher order coefficients TODO: Fix this at some point as well
    for k in 2:order
        # For linear ODEs: x^(k) = J * x^(k-1) / k
        coeffs[k] .= (cache.J * coeffs[k-1]) ./ k
    end
    cache.coeffs = coeffs
    return coeffs
end


function num_signature_matrix(integrator, max_order::Int)
    @unpack u, p, t, f = integrator
    n = length(u)
    Σ = fill(-Inf, n, n)

    residual = zeros(n)
    residual_perturbed = zeros(n)
    du = zero(u)

    compute_residual!(residual, u, du, integrator)
    println("Residual1: ", residual)
    perturbation = 1.1e-8
    for j in 1:n
        for order in 0:max_order
            u_perturbed = copy(u)
            du_perturbed = copy(du)

            if order == 0
                u_perturbed[j] += perturbation
            elseif order == 1
                du_perturbed[j] += perturbation
            else
                continue
            end
            println("Residual: ", residual)
            compute_residual!(residual_perturbed, u_perturbed, du_perturbed, integrator)
            println("Residual perturbed: ", residual_perturbed)
            
            for i in 1:n
                sensitivity = abs(residual_perturbed[i] - residual[i])
                if sensitivity > 1e-8
                    Σ[i, j] = max(Σ[i, j], order)
                end
            end
        end
    end
    return Σ
end



function sum_taylor_series!(integrator, cache)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack Σ, c, d, J, xTS, xtrial, htrial, e, tmp = cache

    order = 5 # This should be fixed or stored or something...
    if hasproperty(integrator.alg, :order)
        order = integrator.alg.order
    end
    coeffs = cache.coeffs
    n = length(u)
    xTS .= uprev 
    factorial_term = 1.0
    h_power = 1.0
    
    for k in eachindex(coeffs)
        h_power *= htrial  # h^k
        factorial_term = k == 1 ? 1.0 : factorial_term * k  # k!
        
        # Add the k-th term: x^(k) * h^k / k!
        xTS .+= coeffs[k] .* (h_power / factorial_term)
    end
    
    # Error estimate
    ets = 0.0
    if length(coeffs) > 0
        # Use the last coefficient to estimate the error
        ets = norm(coeffs[end]) * (htrial^(length(coeffs)+1)) / factorial(length(coeffs)+1)
    end
    
    return xTS, ets
end

using LinearAlgebra
function project!(xtrial, xTS, J, integrator, cache)
    xtrial .= xTS
    
    if integrator.f isa ODEFunction
        # No projection needed for ODEs
        return xtrial
    end
    
    residual = similar(xtrial)
    du = zero(xtrial) # or use cache.du if available

    # Newton iterations to enforce constraints
    for iter in 1:5
        compute_residual!(residual, xtrial, du, integrator)
        if norm(residual) < 1e-10
            break
        end
        
        # TODO: use proper linear solve here (e.g., using Jacobian J)
        xtrial .-= residual  # simple update step
    end

    return xtrial
end


function compute_error!(e_out, xtrial, xTS, ets)
    e = ets + norm(xtrial - xTS)  # Combine Taylor series error and projection error
    e_out = e
    # println("Combined error estimate: ", e)
    return e
end

@muladd function perform_step!(integrator, cache::DAETSCache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack Σ, c, d, J, xTS, xtrial, htrial, e, tmp, xcur = cache

    ttrial = t + htrial
    compute_taylor_coefficients!(integrator, cache)
    xTS, ets = sum_taylor_series!(integrator, cache)
    project!(cache.xtrial, cache.xTS, cache.J, integrator, cache)
    err = compute_error!(cache.e, cache.xtrial, cache.xTS, ets)
    tol = 1e-6  # Set some tolerance
    if hasproperty(integrator, :opts) && hasproperty(integrator.opts, :reltol)
        tol = integrator.opts.reltol
    end
    
    # println("Using tolerance: ", tol)
    # println("Current error: ", err)
    if err <= tol
        # println("Step accepted")
        integrator.u .= cache.xtrial
        integrator.t = ttrial
        if hasproperty(integrator, :tprev)
            integrator.tprev = t
        end
        if hasproperty(integrator, :tcur)
            integrator.tcur = ttrial
        end
        if hasproperty(integrator, :uprev)
            integrator.uprev .= u
        end
        integrator.dt = cache.htrial
    else
        println("Step rejected, adjusting step size")
        # Temporary step size adjust
        new_htrial = htrial * 0.5 
        cache.htrial = max(new_htrial, 1e-10) 
    end
    if hasproperty(integrator, :stats)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 3)
    end

    return nothing
end
"""
    compute_residual!(residual, u, du, integrator)

Compute the residual of the system at given (u, du, t).

- For ODEs, residual = du - f(u, p, t).
- For DAEs, residual is F(residual, du, u, p, t).
"""
function compute_residual!(residual, u, du, integrator)
    @unpack f, p, t = integrator
    if integrator.f isa ODEFunction
        # ODE residual: du - f(u,p,t) = 0
        tmp = similar(u)
        f(tmp, u, p, t)
        @. residual = du - tmp
    else
        # DAE residual: f(residual, du, u, p, t)
        f(residual, du, u, p, t)
    end
    return residual
end
