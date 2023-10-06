using Plots

#variables and constants
l=0.3#metre
g_over_l=5.72

# Defining dimensionless parameters
k_over_msqrtl_over_g = 0.5  # k/m(sqrt(l/g)) 
ω0l_over_g = 2/3  # ω0/g ratio
τ = 3π/100 # Time interval
t_max = 1000 * τ # Total elapsed time

# Initial conditions in dimensionless units
θ0 = 0.0  # Initial angle (θ)
ω0bar = 2/3 # # Initial angular velocity (ω)=dθ/dt


ωbar=2.0 #approximate value


# Defining the function to calculate the angular acceleration (α)
function angular_acceleration(θ::Float64, ωbar::Float64,tbar::Float64, fd0_over_mg::Float64,k_over_msqrtl_over_g::Float64,ω0l_over_g::Float64)::Float64
    α = - sin(θ) - k_over_msqrtl_over_g*ωbar- fd0_over_mg * cos(ω0l_over_g*tbar)
    return α
end

# Implementing the fourth-order Runge-Kutta method to solve the differential equation
function runge_kutta_solver(θ0::Float64, ω0bar::Float64, t_max::Float64, τ::Float64, k_over_msqrtl_over_g::Float64, ω0l_over_g::Float64,fd0_over_mg::Float64)::Tuple
    θ = θ0
    ωbar = ω0bar
    time = 0.0:τ:t_max
    θ_values = [θ]
    ω_values = [ωbar]

    for tbar in time[1:end-1]
        k1_θ = τ * ωbar
        k1_ω = τ * angular_acceleration(θ, ωbar ,tbar, fd0_over_mg,k_over_msqrtl_over_g,ω0l_over_g)

        k2_θ = τ * (ωbar + 0.5 * k1_ω)
        k2_ω = τ * angular_acceleration(θ + 0.5 * k1_θ, ωbar + 0.5 * k1_ω, tbar + 0.5 * τ, fd0_over_mg, k_over_msqrtl_over_g,ω0l_over_g)

        k3_θ = τ * (ωbar + 0.5 * k2_ω)
        k3_ω = τ * angular_acceleration(θ + 0.5 * k2_θ, ωbar + 0.5 * k2_ω, tbar + 0.5 * τ, fd0_over_mg,k_over_msqrtl_over_g,ω0l_over_g)

        k4_θ = τ * (ωbar + k3_ω)
        k4_ω = τ * angular_acceleration(θ + k3_θ, ωbar + k3_ω, tbar + τ, fd0_over_mg, k_over_msqrtl_over_g,ω0l_over_g)

        θ += (k1_θ + 2 * k2_θ + 2 * k3_θ + k4_θ) / 6
        ωbar += (k1_ω + 2 * k2_ω + 2 * k3_ω + k4_ω) / 6

        push!(θ_values, θ)
        push!(ω_values, ωbar)
    end

    return time, θ_values, ω_values
end

# Setting the parameters for the simulation
fd0_over_mg=0.0# Amplitude of the driving force
# Solving the differential equation
@time time_rk, θ_values, ω_values = runge_kutta_solver(θ0, ω0bar, t_max, τ, k_over_msqrtl_over_g, ω0l_over_g ,fd0_over_mg)
# Plotting angular velocity (ω) vs. angle (θ) as time evolves
plot(θ_values,ω_values, label="ω(t) vs θ(t)", xlabel="theta", ylabel="omega")
savefig("omega_vs_theta_plot_0mg.png")
fd0_over_mg=0.89

@time time_rk, θ_values, ω_values = runge_kutta_solver(θ0, ω0bar, t_max, τ, k_over_msqrtl_over_g, ω0l_over_g ,fd0_over_mg)

plot(θ_values,ω_values, label="ω(t) vs θ(t)", xlabel="theta", ylabel="omega")
savefig("omega_vs_theta_plot_0.89mg.png")

fd0_over_mg=1.145
@time time_rk, θ_values, ω_values = runge_kutta_solver(θ0, ω0bar, t_max, τ, k_over_msqrtl_over_g, ω0l_over_g ,fd0_over_mg)
plot(θ_values,ω_values, label="ω(t) vs θ(t)", xlabel="theta", ylabel="omega")
savefig("omega_vs_theta_plot_1.145mg.png")