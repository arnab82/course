using LinearAlgebra
using Plots
using SparseArrays
using IterativeSolvers
# Inputs for creating the grid
Nx = 50 
Ny = 50
Lx = 1 
Ly = 1

# Define the Grid
x = range(0, stop=Lx, length=Nx)
y = range(0, stop=Ly, length=Ny)
dx = x[2] - x[1]
dy = y[2] - y[1]
"""
    relaxation_method_pde_elliptical(Nx::Int, Ny::Int, x, y)

Initialize matrices and boundary conditions for the elliptical partial differential equation.

# Arguments
- `Nx::Int`: The number of grid points along the x-axis.
- `Ny::Int`: The number of grid points along the y-axis.
- `x::Vector{Float64}`: A vector containing the x-coordinate values.
- `y::Vector{Float64}`: A vector containing the y-coordinate values.

# Returns
- `M::Matrix{Float64}`: The coefficient matrix representing the discretized elliptical PDE.
- `B::Vector{Float64}`: The right-hand side vector representing boundary conditions.
"""

function relaxation_method_pde_elliptical(Nx::Int, Ny::Int,x::Vector{Float64},y::Vector{Float64})
        # Initialize matrices
    N = Nx * Ny
    M = zeros(Float64, N, N)
    B = zeros(Float64, N)
    # Interior points
    for i in 2:Nx-1
        for j in 2:Ny-1
            n = i + (j - 1) * Nx
            M[n, n] = -4
            M[n, n - 1] = 1
            M[n, n + 1] = 1
            M[n, n - Nx] = 1
            M[n, n + Nx] = 1
            B[n] = 0
        end
    end

    # Boundary conditions
    # Top BC φ(x, y=∞) = 0, 0 ≤ x ≤ 1
    i = 1
    for j in 1:Ny
        n = i + (j - 1) * Nx
        M[n, n] = 1
        B[n] = 0
    end

    # Bottom BC φ(x=0, y) = 0
    i = Nx
    for j in 1:Ny
        n = i + (j - 1) * Nx
        M[n, n] = 1
        B[n] = 0
    end

    # Left BC φ(x, y=0) = x(1 − x)
    j = 1
    for i in 1:Nx
        n = i + (j - 1) * Nx
        M[n, n] = 1
        B[n] = x[i] * (1 - x[i])
    end

    # Right BC φ(x=1, y) = 0, 0 ≤ y ≤ ∞
    j = Ny
    for i in 1:Nx
        n = i + (j - 1) * Nx
        M[n, n] = 1
        B[n] = 0
    end
    return M,B
end
M,B=relaxation_method_pde_elliptical(Nx,Ny,x,y)
#using the library to see the custom lu decomposition is working fine
# M_sparse = sparse(M)
# F = lu(M_sparse)
# phi_vec = F \ B

#custom lu decomposition
"""
    lu_decomposition_pivoting(A::Matrix)

Perform LU decomposition of a matrix A with partial pivoting.

# Arguments
- `A::Matrix{Float64}`: The input matrix to be decomposed.

# Returns
- `L::Matrix{Float64}`: The lower triangular matrix.
- `U::Matrix{Float64}`: The upper triangular matrix.
- `P::Matrix{Float64}`: The permutation matrix.
"""
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
# Solve for the potential points M*phi=B
@time L,U,P,phi_vec= lu_decomposition(M,lu_decomposition_pivoting,B)

display(phi_vec)

phi = zeros(Float64, Nx, Ny)
for i in 1:Nx
    for j in 1:Ny
        n = i + (j - 1) * Nx
        phi[i, j] = phi_vec[n]
    end
end
display(phi)

heatmap(x, y, phi', c=:viridis, xlabel="x", ylabel="y", aspect_ratio=1)
plot!
savefig("potential_plot.png")