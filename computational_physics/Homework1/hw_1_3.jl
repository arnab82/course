using Printf
function f(x::Float64)::Float64
    return exp(-x)
end

function analytical_integral(a::Float64, b::Float64)::Float64
    return exp(-a) - exp(-b)
end

# Simpson's rule for numerical integration
function simpsons_rule(f::Function, a::Float64, b::Float64, n::Int64)::Float64
    h = (b - a) / n
    x_values=collect(range(a, b, length=n+1))
    # x_values = LinRange(a, b, n + 1)
    integral = 0.0
    for i in 1:n ÷ 2
        integral += (h / 3) * (f(x_values[2*i - 1]) + 4 * f(x_values[2*i]) + f(x_values[2*i + 1]))
    end
    if n % 2 == 0
        integral+=0.0
    else
        integral += (h / 12) * (-f(x_values[n- 1]) + 8 * f(x_values[n]) + 5 * f(x_values[n + 1]))
    end
    
    return integral
end
# a and b are the integration limits
a = 0.0
b = 1.0
exact_integral = analytical_integral(a, b)
#randomly generating the n so that it verifies that the function works for both odd and even numbers
n_list=[998,999,1000,1001,10000,9999]
random_index = rand(1:length(n_list))
n = n_list[random_index]
println("Randomly generated n :",n)
numerical_integral = simpsons_rule(f, a, b, n)
numerical_uncertainty= abs(exact_integral - numerical_integral)

println("Analytical Integral: $exact_integral")
println("Numerical Integral (Simpson's Rule): $numerical_integral")
println("Numerical uncertainty: $numerical_uncertainty")
@printf(" n = %4.2f Analytical Integral =%16.12f numerical integral =%16.12f numerical uncertainty=%16.15f  \n", 10,analytical_integral(a, b),simpsons_rule(f, a, b, 10),abs(analytical_integral(a,b) - simpsons_rule(f, a, b, 10)))
@printf(" n = %4.2f Analytical Integral =%16.12f numerical integral =%16.12f numerical uncertainty=%16.15f  \n", 100,analytical_integral(a, b),simpsons_rule(f, a, b, 100),abs(analytical_integral(a,b) - simpsons_rule(f, a, b, 100)))
@printf(" n = %4.2f Analytical Integral =%16.12f numerical integral =%16.12f numerical uncertainty=%16.15f  \n", 1000,analytical_integral(a, b),simpsons_rule(f, a, b, 1000),abs(analytical_integral(a,b) - simpsons_rule(f, a, b, 1000)))
@printf(" n = %4.2f Analytical Integral =%16.12f numerical integral =%16.12f numerical uncertainty=%16.15f  \n", 10000,analytical_integral(a, b),simpsons_rule(f, a, b, 10000),abs(analytical_integral(a,b) - simpsons_rule(f, a, b, 10000)))
