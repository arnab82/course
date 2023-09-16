using Printf
#a

# Function for nth order polynomial interpolation
function nth_order_interpolation(x::Float64, x_data::Array, y_data::Array)::Float64
    n = length(x_data)-1
    interpolated_value = 0.0
    # f(x)=\sum_{k}^{n+1}f_k*P_k^{n}(x) 
    for i in 1:n+1
        term = y_data[i]
        for j in 1:n+1
            if i != j
                term *= (x - x_data[j]) / (x_data[i] - x_data[j])
            end
        end
        interpolated_value += term
    end
    
    return interpolated_value
end
function error_for_nth_polynomial(x::Float64, x_data::Array, y_data::Array)::Float64
    n = length(x_data)-1
    error=0.0
    #\delta f(x)=f^{n+1}(x)*(x-x1)(x-x2)..............(x-x_{n+1})/(n+1)!
    term=nth_order_interpolation(x,x_data,y_data)/factorial(n+1)
    for i in 1:n+1
        term*=(x-x_data[i])
    end
    error+=term
    return error
end

x_datas=[0.2, 0.3, 0.4, 0.5, 0.6]
y_datas = [0.0015681, 0.0077382, 0.023579, 0.054849, 0.10696]

@printf(" polynomial value = %16.12f error :%16.12f  \n", nth_order_interpolation(0.46,x_datas,y_datas),error_for_nth_polynomial(0.46,x_datas,y_datas))
@printf(" polynomial value = %16.12f error :%16.12f  \n", nth_order_interpolation(0.16,x_datas,y_datas),error_for_nth_polynomial(0.16,x_datas,y_datas))

#b 

function neville_lagrange_interpolation(x::Float64, xdata::Array, ydata::Array)
    n = length(xdata)-1
    function_value= zeros(Float64, n+1, n+1)
    function_value[:, 1] .= ydata
    #neville's interpolation P(x) = \frac{(x - x_j) P_{0,1,\cdots,j-1,j+1,\cdots,k}(x) - (x - x_i) P_{0,1,\cdots,i-1,i+1,\cdots,k}(x)}{(x_i - x_j)}
    for i in 2:n+1
        for j in i:n+1
            function_value[j, i] = ((x - xdata[j - i + 1]) * function_value[j, i - 1] - (x - xdata[j]) * function_value[j - 1, i - 1]) / (xdata[j] - xdata[j - i + 1])
        end
    end
    
    approximated_value = function_value[n+1, n+1]
    
    return Dict("Approximated value" => approximated_value, "Neville iterations table" => function_value)
end

x1=0.16
x2=0.46
result = neville_lagrange_interpolation(x1, x_datas,y_datas)
println("Approximated value at x = $x1: ", result["Approximated value"])
println("Neville iterations table:")
display(result["Neville iterations table"])
result_2= neville_lagrange_interpolation(x2,x_datas,y_datas)
println("Approximated value at x = $x2: ", result_2["Approximated value"])
println("Neville iterations table:")
display(result_2["Neville iterations table"])

function error_neville_lagrange_interpolation(x::Float64, x_data::Array, y_data::Array)::Float64
    result=neville_lagrange_interpolation(x,x_data,y_data)
    n=length(x_data)-1
    numerical_uncertainty=0.5*(abs(result["Neville iterations table"][n+1,n+1]-result["Neville iterations table"][n,n])
                            +abs(result["Neville iterations table"][n+1,n+1]-result["Neville iterations table"][n+1,n]))
    return numerical_uncertainty
end
@printf(" polynomial value =%16.12f in neville's method at x  = %4.2f numerical uncertainty :%16.12f  \n", result["Approximated value"],x1,error_neville_lagrange_interpolation(x1,x_datas,y_datas))
@printf(" polynomial value =%16.12f in neville's method at x  = %4.2f numerical uncertainty :%16.12f  \n", result_2["Approximated value"],x2,error_neville_lagrange_interpolation(x2,x_datas,y_datas))

#c

function linear_interpolation(x_values::Array, y_values::Array, x::Float64)::Float64
    # Number of data points
    n = length(x_values)
    # Finding  the interval which contains x
    i = 1
    while i <= n
        if x < x_values[i]
            break
        end
        i += 1
    end
    # If x is outside the data range, return 0.0
    if i == 1 || i == n + 1
        println("x is outside the data range and hence the function will show 0 as interpolated value")
        return 0.0
    end
    # Linear interpolation
    y = y_values[i - 1] + (y_values[i] - y_values[i - 1]) * (x - x_values[i - 1]) / (x_values[i] - x_values[i - 1])
    return y
end
@printf(" Intrapolated value =%16.12f in linear_interpolation method at x  = %4.2f  \n", linear_interpolation(x_datas,y_datas,x1),x1)
@printf(" Intrapolated  value =%16.12f in linear_interpolation method at x  = %4.2f   \n", linear_interpolation(x_datas,y_datas,x2),x2)
