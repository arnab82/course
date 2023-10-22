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

Nx, Ny = 50, 50
phi = zeros(Float64, Nx, Ny)
h = 0.1

# Values of relaxation parameter ω
omegas = [1.0, 1.2, 1.4, 1.6]
"""
    relaxation_method_analysis(omegas, max_iterations, Nx, Ny, phi, h)

Perform a relaxation method analysis with different relaxation parameters (ω).

# Arguments
- `omegas::AbstractVector`: A vector of relaxation parameters (ω) to analyze.
- `max_iterations::Int`: The maximum number of iterations to perform for each ω.
- `Nx::Int`: The number of grid points along the x-axis.
- `Ny::Int`: The number of grid points along the y-axis.
- `phi::Matrix{Float64}`: The initial field.
- `h::Float64`: The lattice spacing.

# Returns
- `E_values::Dict`: A dictionary containing E values for different ω.
"""
function relaxation_method_analysis(omegas::AbstractVector, max_iterations::Int, Nx::Int, Ny::Int, phi::Matrix{Float64}, h::Float64)
    E_values = Dict()
    for omega in omegas
        E_history = Float64[]
        for iter in 1:max_iterations
            for i in 2:Nx-1
                for j in 2:Ny-1
                    phi[i, j] = (1 - omega) * phi[i, j] + (omega / 4) * (phi[i+1, j] + phi[i-1, j] + phi[i, j+1] + phi[i, j-1])
                end
            end
            E = calculate_energy(phi, h)
            push!(E_history, E)
        end
        E_values[omega] = E_history
    end
    return E_values
end

max_iterations = 1000
E_values=relaxation_method_analysis(omegas,max_iterations,Nx,Ny,phi,h)
plot(1:max_iterations, [E_values[omega] for omega in omegas], label=omegas, xlabel="Iteration Time", ylabel="E", legend=:topright,dpi=300)

savefig("E_vs_iteration_no.png")
