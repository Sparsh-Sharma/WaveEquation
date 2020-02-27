module WaveEquation

# Type to hold the coefficients representing a discretized function
struct DiscretizedFunction{T}
    coeffs::Vector{T}
end

# Arithmetic operations for discretized functions
Base.:+(f::DiscretizedFunction{T}, g::DiscretizedFunction{T}) where T = DiscretizedFunction{T}(f.coeffs + g.coeffs)
Base.:-(f::DiscretizedFunction{T}, g::DiscretizedFunction{T}) where T = DiscretizedFunction{T}(f.coeffs - g.coeffs)
Base.:*(s::Real, f::DiscretizedFunction{T}) where T = DiscretizedFunction{T}(s * f.coeffs)
Base.:*(f::DiscretizedFunction{T}, s::Real) where T = DiscretizedFunction{T}(s * f.coeffs)
Base.:/(f::DiscretizedFunction{T}, s::Real) where T = DiscretizedFunction{T}(1/s * f.coeffs)
Base.:^(f::DiscretizedFunction{T}, s::Real) where T = DiscretizedFunction{T}(f.coeffs.^2)
Base.sum(f::DiscretizedFunction) = sum(f.coeffs)

# Sample an analytically defined function
function discretize(f, N)
    h = 1 / (N-1) # grid spacing
    coeffs = [f(h*i) for i in 0:N-1]
    return DiscretizedFunction(coeffs)
end

# PLC Basis Functions
# gives the value of the ith basis function out of N at position x
# x ranges from 0 to 1
function hat(i, N, x)
    @assert 0 <= i <= N
    h = 1 / (N-1) # grid spacing
    xim1 = h * (i-1)
    xi = h * i
    xip1 = h * (i+1)
    if x <= xim1
        return 0
    elseif x <= xi
        return (x - xim1) / h
    elseif x <= xip1
        return -(x - xip1) / h
    else
        return 0
    end
end

# Reconstruct a function from the discretized coefficients by interpolation
function evaluate(f::DiscretizedFunction{T}, x) where T
    N = length(f.coeffs)
    return sum(f.coeffs[i] * hat(i, N, x) for i in 1:N)
end

# Discrete spatial derivative
function ddx(f::DiscretizedFunction{T}) where T
    N = length(f.coeffs)
    h = 1 / (N-1)
    deriv_coeffs = similar(f.coeffs)
    for i in 1:N
        if i == 1
            deriv_coeffs[i] = (f.coeffs[2] - f.coeffs[1]) / h
        elseif i == N
            deriv_coeffs[i] = (f.coeffs[N] - f.coeffs[N-1]) / h
        else
            deriv_coeffs[i] = (f.coeffs[i+1] - f.coeffs[i-1]) / (2*h)
        end
    end
    return DiscretizedFunction{T}(deriv_coeffs)
end

# Discrete spatial second derivative
function d2dx2(f::DiscretizedFunction{T}) where T
    N = length(f.coeffs)
    h = 1 / (N-1)
    deriv_coeffs = similar(f.coeffs)
    for i in 1:N
        if i == 1
            deriv_coeffs[i] = (f.coeffs[3] - 2*f.coeffs[2] + f.coeffs[1]) / h^2
        elseif i == N
            deriv_coeffs[i] = (f.coeffs[N] - 2*f.coeffs[N-1] + f.coeffs[N-2]) / h^2
        else
            deriv_coeffs[i] = (f.coeffs[i+1] - 2*f.coeffs[i] + f.coeffs[i-1]) / h^2
        end
    end
    return DiscretizedFunction{T}(deriv_coeffs)
end

# Quadrature integration over the domain of f
function quad_integrate(f::DiscretizedFunction{T}) where T
    N = length(f.coeffs)
    h = 1 / (N-1)
    s = T(0)
    for i in 0:N-1
        # if i is at either end, then weight is 0.5, otherwise 1
        w = i==0 || i==N-1 ? T(0.5) : T(1)
        s += w * h * f.coeffs[i+1]
    end
    return s
end

# Runge-Kutta Integration
"""
Given differential equation:
    y'(x) = f(y(x))

Take one RK2 step with:
    f: right hand side function
    y0: state vector y(0)
    h: step size
"""
function rk2step(f, state, dt)
    return state + f(state + f(state)*dt/2)*dt
end

function solve(f, initialstate, dt, n)
    result = [initialstate]
    state = initialstate
    for i in 2:n+1
        newstate = rk2step(f, state, dt)
        push!(result, newstate)
        state = newstate
    end
    return result
end

function waveEQ(state)
    A, Adot = state
    Adot.coeffs[1] = 0
    Adot.coeffs[end] = 0
    return [Adot, d2dx2(A)]
end

export DiscretizedFunction, discretize, hat, evaluate, ddx, d2dx2, rk2step, solve, waveEQ, quad_integrate

end # module
