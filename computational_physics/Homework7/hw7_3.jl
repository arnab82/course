using LinearAlgebra
using Plots

h = 0.05
k = 0.005
x = 0:h:2
t = 0:k:0.05
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

A = diagm(0 => [2 + 2 * factor for i in 1:n-2], 1 => [-factor for i in 2:n-2], -1 => [-factor for i in 2:n-2])
B = diagm(0 => [2 - 2 * factor for i in 1:n-2], 1 => [factor for i in 2:n-2], -1 => [factor for i in 2:n-2])

for j in 1:m-1
    b = T[2:n-1, j] + B * T[2:n-1, j]
    b[1] += factor * (T[1, j] + T[1, j+1])
    b[end] += factor * (T[end, j] + T[end, j+1])
    solution = A \ b
    T[2:n-1, j+1] = solution
end



labels = ["Point $i" for i in t]
# display(T)
plot(x, T, color = :auto,label=labels,linewidth=3,
     xlabel = "distance [m]", ylabel = "Temperature [°C]", legend = false)
savefig("crank_nicolson_diffusion.png") 
