import numpy as np
import matplotlib.pyplot as plt
import time
start_time = time.time()
# Defining  the ODE 
def ode_system(x:float, y:np.ndarray):
    y_, d2ydx2 = y[1], -4 * np.pi**2 * y[0]
    return [y_, d2ydx2]


# Defining the 4th-order Runge-Kutta method
def runge_kutta_step_4th_order(x:float, y:np.ndarray, h:float):
    k1 = np.multiply(h, ode_system(x, y))
    k2 = np.multiply(h, ode_system(x + h / 2, y + k1 / 2))
    k3 = np.multiply(h, ode_system(x + h / 2, y + k2 / 2))
    k4 = np.multiply(h, ode_system(x + h, y + k3))
    # print(np.add(y, (k1 + 2 * k2 + 2 * k3 + k4) / 6))
    return np.add(y, (k1 + 2 * k2 + 2 * k3 + k4) / 6)
    


# Defining the secant method
def custom_secant_method(equation, x0:float, x1:float, tol:float, max_iter:int):
    for i in range(max_iter):
        f0 = equation(x0)
        f1 = equation(x1)
        if abs(f0-f1) < tol:
            return x1  # when Converged to a solution within tolerance
        x_2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        x0, x1 = x1, x_2
    return None  # when it do not converge within maximum no of iterations



# Defining the shooting method with runge-kutta method and secant method as the root finding algorithm
def shooting_method_with_runge_kutta_and_secant(y0_guess:float, dydx0_guess:float, x_range:np.ndarray, boundary_condition_dydx:float):
    def objective(d2ydx20):#defining the objective function
        y0 = [y0_guess, d2ydx20]#initial values of y and its second derivative
        x_values, y_values = [], []
        # Performing the runge-kutta step 
        for x in x_range:
            x_values.append(x)
            y_values.append(y0[0])
            y0 = runge_kutta_step_4th_order(x, y0, x_range[1] - x_range[0])
        # Return the value of dy/dx at x=1
        return y0[1] - boundary_condition_dydx #applying the boundary condition
    # Initial guesses for the secant method
    x0 = dydx0_guess - 0.1
    x1 = dydx0_guess + 0.1
    dydx0_solution = custom_secant_method(objective, x0, x1, tol=1e-6, max_iter=100)
    return dydx0_solution

# Analytical solution
def analytical_solution(x):
    return np.sin(2 * np.pi * x) + np.cos(2 * np.pi * x)




# Defining the x values for the numerical solution
x_range = np.linspace(0, 1, 1000)#if we increase the number of points from 100 to 1000 , we can see the the convergence of numerical solution to analytical solution
# Guess for initial dy/dx at x=0
dydx0_guess = 0.0
# Target dy/dx at x=1
dydx_at_x_eq_1 = 2 * np.pi
# Solving the ODE using the shooting method
dydx0_solution = shooting_method_with_runge_kutta_and_secant(1.0, dydx0_guess, x_range, dydx_at_x_eq_1 )
# Calculating the corresponding solution using the 4th-order Runge-Kutta method
y0_solution = [1, dydx0_solution]
y_solution = []
for x in x_range:
    y_solution.append(y0_solution[0])
    y0_solution = runge_kutta_step_4th_order(x, y0_solution, x_range[1] - x_range[0])
# print(y_solution)
# Analytical solution
y_analytical = analytical_solution(x_range)



# Function for plotting
def plot_solution(x_range:np.ndarray, y_solution:np.ndarray, y_analytical:np.ndarray,filename):
    plt.figure(figsize=(12, 8))
    plt.plot(x_range, y_solution, label="Numerical Solution")
    plt.plot(x_range, y_analytical, label="Analytical Solution", linestyle="--")
    plt.xlabel("x")
    plt.ylabel("y(x)")
    plt.legend()
    plt.title("Shooting Method vs Analytical Solution")
    plt.grid(True)
    plt.savefig(filename)
    plt.show()
plot_solution(x_range,y_solution,y_analytical,"shooting_vs_analytical")
end_time = time.time()
elapsed_time = end_time - start_time
print(f"Execution Time: {elapsed_time} seconds")
