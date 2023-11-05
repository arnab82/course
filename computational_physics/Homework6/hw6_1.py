import numpy as np
from scipy.integrate import odeint
from scipy.integrate import quad
from scipy.optimize import bisect

# Define the differential equation
def diff_eq(y, x, alpha):
    return [y[1], -4 * np.pi**2 * y[0] - alpha * y[0]]

# Step (ii): Initial conditions
x0 = 0.0
u0 = 1.0
h = 1e-3  # Step size

# Calculate u1 using the second point
x1 = 1
u1 = u0 + h * 2 * np.pi  # dy/dx(x=0) = 2π

# Step (iv): Find the root of F(λ)
def F(lambda_val):
    y0 = [u0, lambda_val]  # Initial conditions at x=0
    x_values = np.linspace(x0, x1, 100)
    solution = odeint(diff_eq, y0, x_values, args=(lambda_val,))
    u_lambda = solution[-1, 0]  # Value of u_λ at x=1

    return u_lambda - u1


lower_bound = -10.0 
upper_bound = 10.0
lambda_root = bisect(F, lower_bound, upper_bound)

# Print the results
print(f"Root λ: {lambda_root}")

# Step (v): Calculate the corresponding eigenvector with normalization
def eigenvector(u, integral):
    return u / np.sqrt(integral)

# Calculate the integral using quad
def u_lambda(x, lambda_val):
    y0 = [u0, lambda_val]  # Initial conditions at x=0
    x_values = np.linspace(x0, x1, 100)
    solution = odeint(diff_eq, y0, x_values, args=(lambda_val,))
    integral, _ = quad(lambda x: solution[:, 0]**2, x0, x1)

    return integral
# Calculate the integral using a loop
integral_value = 0.0
num_intervals = 100  # The number of intervals for numerical integration
x_values_integration = np.linspace(0, 1, num_intervals)

for i in range(len(x_values_integration) - 1):
    a = x_values_integration[i]
    b = x_values_integration[i + 1]
    integral, _ = quad(lambda x: u_lambda(x, lambda_root)**2, a, b)
    integral_value += integral

# Calculate the eigenvector with normalization
eigenvector_result = eigenvector(u_lambda(x1, lambda_root), integral_value)


