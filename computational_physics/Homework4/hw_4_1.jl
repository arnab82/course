using LinearAlgebra
import Arpack: eigs
# Define the matrix A and vector b
A = [2.0 3.0 10.0 -1.0;
        10 15 3 7;
        -4 1 2 9;
        0 0 0 0]
b =[1, 2, 3,0]
using LinearAlgebra
function custom_svd(A, tol=1e-10)
    # Step 1: Compute the covariance matrix
    ATA = A' * A
    # println("ATA")
    # display(ATA)
    
    # Step 2: Compute eigenvalues and eigenvectors of ATA
    eigenvalues, eigenvectors = eigen(ATA)
    
    # Sort eigenvalues and eigenvectors in descending order
    sorted_indices = sortperm(eigenvalues, rev=true)
    eigenvalues = eigenvalues[sorted_indices]
    eigenvectors = eigenvectors[:, sorted_indices]
    # println("Eigenvalues")
    # display(eigenvalues)
    # println("Eigenvectors")
    # display(eigenvectors)
    
    # Step 3: Calculate the singular values and sort them
    singular_values = sqrt.(eigenvalues)
    sorted_singular_indices = sortperm(singular_values, rev=true)
    singular_values = singular_values[sorted_singular_indices]
    eigenvectors = eigenvectors[:, sorted_singular_indices]
    # println("Singular values")
    # display(singular_values)
    
    # Step 4: Calculate the matrix V
    V = eigenvectors
    
    U = zeros(size(A, 2), size(A, 2))
    
    # Step 5: Calculate the matrix U
    for i in 1:length(eigenvalues)
        U[:, i] = A * V[:, i] / norm(A * V[:, i])
        # println(U[:,i])
    end
    # Step 6: Ensure that U, V are real
    U = real(U)
    V = real(V)
    
    # Step 7: Handle the case when singular values are close to zero
    non_zero_singular_values = singular_values .> tol
    k = sum(non_zero_singular_values)

    singular_values = singular_values[1:k]
    w = zeros(k, k)
    
    for i in 1:length(singular_values)
        w[i, i] = singular_values[i]
    end
    
    return U, w, transpose(V)
end

# Perform SVD using the custom_svd function
U, w, VT = custom_svd(A)
singular_value_matrix = Diagonal(w)

# Reconstruct A using U, singular_values, and VT
reconstructed_A = U * (singular_value_matrix * VT)
println("U")
display(U)
println("W")
display(w)
println("Singular values")
println(singular_value_matrix)
println("VT")
display(VT)
println("Reconstructed A")
display(reconstructed_A)

# Check if the reconstruction is close to the original matrix A
reconstruction_error = norm(A - reconstructed_A)
println("Reconstruction Error:", reconstruction_error)

# Step 3: Solve for x
# Calculate the pseudo-inverse of W (W^+)
tolerance = 1e-6 # A small value to handle singular values close to zero
w_inv = (ifelse.(singular_value_matrix .> tolerance, 1.0 ./ singular_value_matrix, 0.0))
println("Inverse of W")
display(w_inv)

# Calculate x using SVD components
x = VT' * w_inv * transpose(U) * b

# Step 4: General Solution
# The homogeneous solution (x_h) depends on the null space of A
# We can find the null space using the SVD components

# Calculate the homogeneous solution (x_h)
x_h = VT[:, end]  # Take the last column of VT (corresponding to the smallest singular value)

# Step 5: Confirm the Solution
# Check if Ax_p is close to b
residual = norm(A * x - b)
println("Residual (Ax_p - b):", residual)

# Print the particular solution (x_p)
println("Particular Solution (x_p):", x)

# Print the homogeneous solution (x_h)
println("Homogeneous Solution (x_h):", x_h)


#using linear algebra in built svd function

F= svd(A)
# Display the results
println("U:")
display(F.U)
println("\nS:")
display(F.S)
println("\nVt:")
display(F.Vt)
reconstructed_A_ = F.U * Diagonal(F.S)* F.Vt
println("Reconstructed A")
display(reconstructed_A_)


tolerance = 1e-6  
S_inv = Diagonal([s > tolerance ? 1.0 / s : 0.0 for s in F.S])

x = F.Vt' * S_inv * F.U' * b
x_h = F.Vt[:, end] 
residual = norm(A * x - b)
println("Residual (Ax_p - b): $residual")
println("Particular Solution (x_p):")
println(x)
println("Homogeneous Solution (x_h):")
println(x_h)