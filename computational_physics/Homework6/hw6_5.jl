using LinearAlgebra
using Plots
"""
    calculate_energy(phi, h)

Calculate the energy functional E based on a numerical solution of Laplace's equation.

# Arguments
- `phi::Matrix{Float64}`: A 2D array representing the numerical solution of Laplace's equation.
- `h::Float64`: The lattice spacing used in the discretization.

# Returns
- `E::Float64`: The computed energy based on the discretized energy functional.

"""
function calculate_energy(phi::Matrix{Float64}, h::Float64)
    # Calculate the energy E based on the discretized energy functional
    E = 0.0
    Nx, Ny = size(phi)
    for i in 2:Nx-1
        for j in 2:Ny-1
            gradient_term = (phi[i+1, j] + phi[i-1, j] + phi[i, j+1] + phi[i, j-1] - 4 * phi[i, j]) / h^2
            E += gradient_term^2
        end
    end
    E /= h^2  # Normalize by the lattice size
    return E
end
"""
    convergence_analysis(h_values, max_iterations)

Perform a convergence analysis of the energy functional E for different lattice spacings and iterations.

# Arguments
- `h_values::Vector{Float64}`: An array of lattice spacings to analyze.
- `max_iterations::Int`: The maximum number of iterations for each lattice spacing.

# Returns
- `convergence_data::Dict`: A dictionary containing convergence data for each lattice spacing.
"""

function convergence_analysis(h_values::Vector{Float64}, max_iterations::Int)
    convergence_data = Dict()

    for h in h_values
        Nx, Ny = 10, 10
        phi = zeros(Float64, Nx, Ny)

        E_history = Float64[]
        for iter in 1:max_iterations
            # Update phi using relaxation method
            for i in 2:Nx-1
                for j in 2:Ny-1
                    phi[i, j] = 0.25 * (phi[i+1, j] + phi[i-1, j] + phi[i, j+1] + phi[i, j-1])
                end
            end
            E = calculate_energy(phi, h)
            push!(E_history, E)
        end

        convergence_data[h] = E_history
    end

    return convergence_data
end

# Define lattice spacings (h values) for coarse and fine meshes
h_coarse = 0.2
h_fine = 0.05

max_iterations = 1000
h_values = [h_coarse, h_fine]
convergence_data = convergence_analysis(h_values, max_iterations)
plot(1:max_iterations, [convergence_data[h] for h in h_values], label=h_values, xlabel="Mesh size", ylabel="E", legend=:topright,dpi=300)
savefig("E_vs_scalability_of_meshes.png")