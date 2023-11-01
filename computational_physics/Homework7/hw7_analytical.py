import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(0, 2 , 50)
t = 0.045
# Calculate Ψ(x, t) using the provided equation
def psi(x, t):
    result = 0.0
    for n in range(0, len(x)): 
        result += (-1)**n / ((2 * n + 1) ** 2 * (np.pi ** 2)) * np.exp(
                -((2 * n + 1) ** 2 * (np.pi ** 2) * t) / 4) * np.sin((n + 1 / 2) * np.pi * x)
    return result * 8.0
# Calculate Ψ(x, t) for the given x and t
psi_values = psi(x, t)
plt.plot(x, psi_values)
plt.xlabel('x')
plt.ylabel('Ψ(x, t)')
plt.title('Wave Function Ψ(x, t)')
plt.grid(True)
plt.savefig("analytical_solution.png")
plt.show()
