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

# Explicit method
for j in 2:m
    for i in 2:n-1
        T[i, j] = factor * T[i - 1, j - 1] + (1 - 2 * factor) * T[i, j - 1] + factor * T[i + 1, j - 1]
    end
end



matrix1 = []

push!(matrix1, T[:, round(Int, 0.005 / k)])
push!(matrix1, T[:, round(Int, 0.01 / k)])
push!(matrix1, T[:, round(Int, 0.02 / k)])
push!(matrix1, T[:, round(Int, 0.03 / k)])
push!(matrix1, T[:, round(Int, 0.045 / k)])

labels = [0.005, 0.01, 0.02, 0.03, 0.045]

plot(x, [matrix1[1],matrix1[2],matrix1[3],matrix1[4],matrix1[5]], color=:auto, xlabel="distance[m]", ylabel="Temperature[° C]",linewidth=3)

savefig("explicit_diffusion_1.png") 