using LinearAlgebra
using Printf
function lu_decomposition_pivoting(A::Matrix)
    n = size(A, 1)
    L = Matrix{Float64}(I, n, n)#initializing Lower traingular matrix as an identity matrix
    U = copy(A)                 #initializing upper traingular matrix as a copy of the original matrix
    P = Matrix{Float64}(I, n, n)#initializing permutation matrix as an identity matrix

    for k in 1:n-1
        pivot_row = argmax(abs.(U[k:end, k])) + k - 1#finding the row index with the maximum absolute value in the current column k for
                                                    # U matrix and adding the k-1 for giving the actual index in the original U matrix as argmax considers k as the first row.

        if pivot_row != k
            # Swapping rows in U and apply the same permutation to P
            U[[k, pivot_row], :] = U[[pivot_row, k], :]
            P[[k, pivot_row], :] = P[[pivot_row, k], :]

            if k > 1
                # Appling the same permutation to L for rows below the first row as it contains only one non-zero element(1) in the diagonal
                L[[k, pivot_row], 1:k-1] = L[[pivot_row, k], 1:k-1]
            end
        end

        for i in k+1:n
            factor = U[i, k] / U[k, k]
            L[i, k] = factor
            U[i, k:end] .-= factor * U[k, k:end]
        end
    end

    return L, U, P
end


function lu_decomposition(A,lu_decomposition_pivoting::Function, b::Vector)
    L,U,P=lu_decomposition_pivoting(A)
    # Solve Ly = Pb
    y = L \ (P * b)

    # Solve Ux = y
    x = U \ y
    return L,U,P,x
end

# Given matrix A and right-hand side vector b
A = [1.0 -1.0 -2.0 1.0;
     3.0 2.0 -1.0 2.0;
     2.0 3.0 1.0 -3.0;
     10.0 -4.0 3.0 2.0]

b = [1.0, 4.0, 2.0, 3.0]
# Getting the L,U,P and solution vector
@time L,U,P,x= lu_decomposition(A,lu_decomposition_pivoting,b)
# Display L, U, and P matrices
println("L Matrix:")
display(L)
println("U Matrix:")
display(U)
println("P Matrix:")
display(P)
# Checking if L·U = P·A
PA = P * A
LU = L * U

println("LU Matrix:")
display(LU)
println("PA Matrix:")
display(PA)
println("Checking if L·U = P·A")
println("L·U == P·A:\n", isapprox(LU, PA))

# Display the solution vector
println("Solution vector")
@printf(" x\t y\t  z\t u\n%6.4f %6.4f %6.4f %6.4f\n", x[1], x[2], x[3], x[4])

println("Check if the numerical solution satisfies the given answer")
given_answer = A \ b
println("Numerical solution matches the given answer:", isapprox(x, given_answer))

# Computing the determinant of matrix A
det_A = det(A)
println("Determinant of matrix A:", det_A)
