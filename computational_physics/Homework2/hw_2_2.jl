function gauss_elimination(A::Matrix, b::Vector)::Vector
    n = length(b)
    Ab = [A b]# stacking the augmented matrix and the vector
    
    # Forward elimination
    for k in 1:n
        # Finding the pivot row and swapping rows if required ; although pivoting is not required for this particular problem
        maxval, pivotrow = findmax(abs.(Ab[k:n, k]))
        pivotrow += k - 1
        Ab[k, :], Ab[pivotrow, :] = Ab[pivotrow, :], Ab[k, :]
        # display(Ab)
        # Eliminate the entries below the pivot
        for i in k+1:n
            factor = Ab[i, k] / Ab[k, k]
            # println(factor)
            # display(Ab[i,k:n+1])
            # display(Ab[k,k:n+1])
            # display(factor .* Ab[k, k:n+1])
            Ab[i, k:n+1] = Ab[i, k:n+1].-factor .* Ab[k, k:n+1]
        end
    end
    
    # Backward substitution
    x = zeros(Float64, n)
    for i in n:-1:1
        x[i] = (Ab[i, end] - sum(Ab[i, i+1:n].* x[i+1:n])) / Ab[i, i]
    end
    
    return x
end

# Coefficient matrix A 
A = [10.0 -4.0 3.0 2.0; 1.0 -1.0 -2.0 1.0; 3.0 2.0 -1.0 2.0; 2.0 3.0 1.0 -3.0]
display(A)
b = [3.0, 1.0, 4.0, 2.0]
display(b)

# Solving the system of equations
@time x = gauss_elimination(A, b)
println("Solution:")
println("x = ", x[1])
println("y = ", x[2])
println("z = ", x[3])
println("u = ", x[4])
