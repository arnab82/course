using LinearAlgebra
using Polynomials
using Plots
"""
    calculate_energy(h)

Calculate the energy functional E based on the relaxation method for Laplace's equation.

# Arguments
- `h::Float64`: The lattice spacing.

# Returns
- `E::Float64`: The computed energy based on the discretized energy functional.
"""

function calculate_energy(h)
    Nx, Ny = 10, 10  # Define the grid size
    phi = zeros(Float64, Nx, Ny)  # Initialize the phi matrix

    # Relaxation method for updating phi
    for iter in 1:1000 
        # Iterate and update phi using relaxation method
        for i in 2:Nx-1
            for j in 2:Ny-1
                phi[i, j] = 0.25 * (phi[i+1, j] + phi[i-1, j] + phi[i, j+1] + phi[i, j-1])
            end
        end
    end

    # Calculate the energy based on the discretized energy functional
    E = 0.0
    for i in 2:Nx-1
        for j in 2:Ny-1
            gradient_term = (phi[i+1, j] + phi[i-1, j] + phi[i, j+1] + phi[i, j-1] - 4 * phi[i, j]) / h^2
            E += gradient_term^2
        end
    end

    # Normalize by the lattice size
    E /= h^2

    return E
end

# Array of lattice spacings
h_values = [0.1, 0.05, 0.01, 0.005]
E_values = [calculate_energy(h) for h in h_values]

# Fit the data  a quadratic polynomial
degree = 2
fit = polyfit(h_values, E_values, degree)
E_h0 = fit[end]
println("E(h → 0) = $E_h0")
scatter(h_values, E_values, label="Data Points", xlabel="h", ylabel="E")
plot!(h_values, polyval(fit, h_values), label="Polynomial Fit",dpi=300)
savefig("energy_vs_h.png")
