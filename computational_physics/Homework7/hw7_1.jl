using Plots

h = 0.045
k = 0.0015
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

# Explicit method
for j in 2:m
    for i in 2:n-1
        T[i, j] = factor * T[i - 1, j - 1] + (1 - 2 * factor) * T[i, j - 1] + factor * T[i + 1, j - 1]
    end
end
labels = ["Point $i" for i in t]
# display(T)
plot(x, T, color = :auto,label=labels,linewidth=3,
     xlabel = "distance [m]", ylabel = "Temperature [°C]", legend = true)
savefig("explicit_diffusion.png") 