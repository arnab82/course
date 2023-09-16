using Printf


function secant_method(equation::Function, a::Float64, x_values::Array, tol::Float64, max_iter::Int64)
    random_index1 = rand(1:length(x_values))
    random_index2 = rand(1:length(x_values))
    x_i= x_values[random_index1]
    println(x_i)
    x_j=x_values[random_index2]
    println(x_j)
    if equation(x_i,a)*equation(x_j,a)<0.0# initial condition to verify if the points are applicable for secant methods
        for i in 1:max_iter
            f_i = equation(x_i, a)
            f_j= equation(x_j, a)
            if abs(x_j-x_i) > tol 
                x_next = x_j- f_j * (x_j - x_i) / (f_j - f_i)# computing the next guess for x 
                x_i, x_j = x_j, x_next  
                @printf(" Iteration No: %4i Abs(xj-xi): %16.12f Xj: %16.12f\n", i,abs(x_j-x_i),x_j )
            elseif i==max_iter
                println("it didn't converge")
            else
                return x_j# Converged to a solution within tolerance
                break 
            end
        end
    
    else
        println("the initial points are not suitable for this method")    
    end 
end

# Example usage:
equation(x, a) = tan(x) - a / x  # Replace with your equation
a = 5.0
tolerance = 1e-8
max_iterations = 100
x_values=collect(range(pi/2+tolerance,3*pi/2-tolerance, length=20))# limit value of x is not taken as 
                                                                    #tan(pi/2)and tan(3*pi/2) goes to infinity ; so the roots will not be accurate
# creating random initial guess x value in the limit
random_index1 = rand(1:length(x_values))
random_index2 = rand(1:length(x_values))
x_i= x_values[random_index1]
x_j=x_values[random_index2]


result = secant_method(equation, a, x_values, tolerance, max_iterations)
if result !== nothing
    println("Approximated root: $result")
else
end
random_result=[]

@time for i in 1:20
    @time result1 = secant_method(equation, a, x_values, tolerance, max_iterations)
    push!(random_result,result1)
end 
println(" \n\n The list of approximated roots for 20 random initial guess; \nif they are not appropriate points for secant method, nothing will be written in list")
display(random_result)  
