using LinearAlgebra
using Plots

h = 0.04
k = 0.00075
x = 0:h:2
t = 0:k:0.045
boundary_conditions = [0, 0]
n = length(x)
m = length(t)
T = zeros(n, m)
T[1, :] .= boundary_conditions[1]
T[end, :] .= boundary_conditions[2]

# Initial conditions
for i in 1:n
    if 0 <= x[i] <= 1
        T[i, 1] = x[i]
    else
        T[i, 1] = -x[i] + 2
    end
end

factor = k / h^2

A = diagm(0 => [2 + 2 * factor for i in 1:n-2], 1 => [-factor for i in 2:n-3], -1 => [-factor for i in 2:n-3])
B = diagm(0 => [2 - 2 * factor for i in 1:n-2], 1 => [factor for i in 2:n-3], -1 => [factor for i in 2:n-3])

for j in 1:m-1
    b = deepcopy(T[2:n-1, j])
    b= B * b
    b[1] += factor * (T[1, j] + T[1, j+1])
    b[end] += factor * (T[end, j] + T[end, j+1])
    solution = A \ b
    print(solution)
    T[2:n-1, j+1] = solution
end


matrix1 = []

push!(matrix1, T[:, round(Int, 0.005 / k)])
push!(matrix1, T[:, round(Int, 0.01 / k)])
push!(matrix1, T[:, round(Int, 0.02 / k)])
push!(matrix1, T[:, round(Int, 0.03 / k)])
push!(matrix1, T[:, round(Int, 0.045 / k)])

labels = [0.005, 0.01, 0.02, 0.03, 0.045]

plot(x, [matrix1[1],matrix1[2],matrix1[3],matrix1[4],matrix1[5]], color=:auto, xlabel="distance[m]", ylabel="Temperature[° C]",linewidth=3)
savefig("crank_nicolson_diffusion_2.png") 
t1 = 0.045

function psi_analytical(x, t)
    result = 0.0
    for n in 0:length(x)-1
        result += (-1)^n / ((2n + 1)^2 * (π^2)) * exp(-((2n + 1)^2 * (π^2) * t) / 4) * sin((n + 1/2) *22/7* x)
    end
    return result * 8
end

psi_values = psi_analytical.(x, t1)
numerical_values = matrix1[5]

plot(x, psi_values, label="Analytical", color="red")
plot!(x, numerical_values, label="Numerical", color="orange")

xlabel!("x")
ylabel!("Ψ(x, t)")
title!("Wave Function Ψ(x, t) : analytical vs numerical")


savefig("analytical_solution_vs_numerical_crank_nicholson.png")

