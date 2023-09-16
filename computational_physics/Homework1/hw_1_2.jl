using Printf

# Function to calculate f(x)
function f(x::Float64)::Float64
    return x^3 * cos(x)
end

# First-order derivative using three-point formula
function first_derivative(x::Float64, h::Float64)::Float64
    return (f(x + h) - f(x - h)) / (2 * h)
end

# Second-order derivative using three-point formula
function second_derivative(x::Float64, h::Float64)::Float64
    return (f(x + h) - 2 * f(x) + f(x - h)) / (h^2)
end

# Analytical first and second derivatives
function analytical_first_derivative(x::Float64)::Float64
    return 3x^2 * cos(x) - x^3 * sin(x)
end

function analytical_second_derivative(x::Float64)::Float64
    return 6x * cos(x) - 6x^2 * sin(x) - x^3 * cos(x)
end


x_start = π / 2
x_end = π
num_intervals = 100
h = (x_end - x_start) / num_intervals
x_values = Float64[]
numerical_first_derivatives = Float64[]
numerical_second_derivatives = Float64[]
analytical_second_derivatives = Float64[]
analytical_first_derivatives = Float64[]
numerical_errorfirst = Float64[]
numerical_errorsecond = Float64[]
for i in 0:num_intervals
    x = x_start + i * h
    
    # Calculating derivatives using three-point formulas
    first_deriv = first_derivative(x, h)
    second_deriv = second_derivative(x, h)
    
    # Calculating analytical derivatives
    analytical_first = analytical_first_derivative(x)
    analytical_second = analytical_second_derivative(x)
    
    # Calculating numerical accuracy
    numerical_error_first = abs(analytical_first - first_deriv)
    numerical_error_second = abs(analytical_second - second_deriv)

    push!(x_values, x)
    push!(numerical_first_derivatives, first_deriv)
    push!(numerical_second_derivatives, second_deriv)
    push!(numerical_errorfirst, numerical_error_first)
    push!(numerical_errorsecond, numerical_error_second)
end

println("x\t\tNumerical 1st Derivative\tNumerical 2nd Derivative\tNumerical Error (1st)\tNumerical Error (2nd)")
for i in 1:length(x_values)
    @printf("%.6f\t\t%.8f\t\t%.8f\t\t\t\t%.6f\t\t%.6f\n", x_values[i], numerical_first_derivatives[i], numerical_second_derivatives[i], numerical_errorfirst[i], numerical_errorsecond[i])
end
