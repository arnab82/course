using Plots

# Define the domain and discretization
Nx, Ny = 50, 50
Lx, Ly = 1, 1
dx, dy = Lx / (Nx - 1), Ly / (Ny - 1)

# Define the grid
x = range(0.0, stop=Lx, length=Nx)
y = range(0.0, stop=Ly, length=Ny)

# Initialize the phi matrix
phi = zeros(Nx, Ny)

# Set the boundary conditions
for i in 1:Nx
    phi[i, 1] = (i - 1) * dx * (1.0 - (i - 1) * dx)  # φ(x, y = 0) = x(1 - x)
    phi[i, Ny] = 0.0  # φ(x, y = ∞) = 0
end
for j in 1:Ny
    phi[1, j] = 0.0  # φ(x = 0, y) = 0
    phi[Nx, j] = 0.0  # φ(x = 1, y) = 0
end


num_iterations = 1000
# Initialize a plot
heatmap(x, y, phi', c=:viridis, xlabel="x", ylabel="y", aspect_ratio=1)
"""
    perform_jacobi_iteration(phi, num_iterations)

Perform Jacobi iterations to solve a Laplace equation.

# Arguments
- `phi::Matrix{Float64}`: The initial field.
- `num_iterations::Int`: The number of Jacobi iterations to perform.

# Returns
- `phi::Matrix{Float64}`: The updated field after the Iteration
"""
function perform_jacobi_iteration(phi::Matrix{Float64}, num_iterations::Int,x::AbstractVector, y::AbstractVector)
    Nx, Ny = size(phi)
    new_phi = copy(phi)
    for iter in 1:num_iterations
        for i in 2:Nx-1
            for j in 2:Ny-1
                new_phi[i, j] = 0.25 * (phi[i+1, j] + phi[i-1, j] + phi[i, j+1] + phi[i, j-1])
            end
        end
        phi .= new_phi
        # Update the plot for every 1 iterations 
        if iter % 1 == 0
            plot!(x, y, phi', c=:viridis,dpi=300)
        end
    end
    return phi
end
perform_jacobi_iteration(phi,1000,x,y)
display(plot)
savefig("potential_plot_analytical.png")