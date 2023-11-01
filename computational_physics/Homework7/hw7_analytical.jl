using Plots

x = range(0, stop=2, length=10000)
t = 0.045

function psi(x, t)
    result = 0.0
    for n in 0:length(x)-1
        result += (-1)^n / ((2n + 1)^2 * (π^2)) * exp(-((2n + 1)^2 * (π^2) * t) / 4) * sin((n + 1/2) *22/7* x)
    end
    return result * 8
end

psi_values = psi.(x, t)

plot(x, psi_values, xlabel="x", ylabel="Ψ(x, t)", title="Wave Function Ψ(x, t)", legend=false,linewidth=3)

savefig("analytical_solution.png")
