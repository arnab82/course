using Printf
using LinearAlgebra
function generalised_partial_pivoting(A::Matrix, b::Vector)
    n = size(A)[1]
    for col in 1:n
        # Finding the row with the maximum value in the current column(this will search each column but only for lower triangular matrix i.e below the the diagonal elements)
        max_row = argmax(abs.(A[col:end, col])) + col - 1
        
        if max_row != col
            # Swap rows in A, b
            A[[col, max_row], :] = A[[max_row, col], :]
            b[[col, max_row]] = b[[max_row, col]]
        end
    end
    
    return A, b
end


function gauss_elimination(A::Matrix, b::Vector)
    n = size(A)[1]
    
    # Forward elimination
    for i in 1:n#columnwise
        for j in i+1:n#rowwise
            A, b= generalised_partial_pivoting(A, b)
            factor = A[j, i] / A[i, i]
            A[j, i:end] .-= factor * A[i, i:end]
            b[j] -= factor * b[i]
            # display(A)
            display(hcat(A,b))
        end
    end
    
    # Backward substitution
    x = zeros(n)
    for i in n:-1:1#We are moving from the last row to upward direction
        x[i] = (b[i] - dot(A[i, i+1:end], x[i+1:end])) / A[i, i]
    end
    
    return x
end

# Given matrix A and right-hand side vector b
A = [1.0 -1.0 -2.0 1.0;
     3.0 2.0 -1.0 2.0;
     2.0 3.0 1.0 -3.0;
     10.0 -4.0 3.0 2.0]

b = [1.0, 4.0, 2.0, 3.0]

# Solving the system of equations via gauss elimination
@time x = gauss_elimination(A, b)
expected_solution = [21/34, 43/68, -13/34, 1/4]

@printf("Evaluated Solution (x, y, z, u): %6.4f %6.4f %6.4f %6.4f\n", x[1], x[2], x[3], x[4])
@printf("Expected Solution (x, y, z, u): %6.4f %6.4f %6.4f %6.4f\n", expected_solution[1], expected_solution[2], expected_solution[3], expected_solution[4])