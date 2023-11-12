import random
import numpy as np
import matplotlib.pyplot as plt

def W(x):
    return (x**2+2.0)

M = 10000# Number of sample points
N = 0.0  # Initialize normalization constant

# Determine the normalization constant N by integrating W(x) over the range [0, 1]
for _ in range(M):
    x = random.uniform(0, 1)
    N += W(x)
N /= M

def integrand(x):
    return x**2

# Perform the Monte Carlo integration with importance sampling
def monte_carlo_integral(M, N):
    integral_sum = 0.0
    total_weight = 0.0  
    sample_means = []
    sample_points = []
    total_points = []

    for i in range(M + 80000):
        x = random.uniform(0, 1)
        weight = W(x) /N  # Importance sampling weight
        integrand_val = integrand(x)
        total_points.append(integrand_val/weight)
        if i >= 80000:
            sample_points.append(integrand_val)
            integral_sum += integrand_val 
            total_weight += weight
            if i >= 80000 and (i + 1) % 100 == 0:
            # Calculate the sample mean for every 100 points
                sample_mean = integral_sum / total_weight
                sample_means.append(sample_mean)
    # Calculate the integral result
    integral_result = integral_sum / total_weight

    return integral_result, sample_means, sample_points,total_points

# Calculate the integral, sample means, and sample points
integral_result, sample_means, sample_points ,total_points= monte_carlo_integral(M, N)

# Calculate the standard deviation of the sample mean
standard_deviation = np.std(sample_means)
print("Monte Carlo Integral Result:", integral_result)
print("Standard Deviation of Sample Mean:", standard_deviation)


# \begin{array}{l}
# C(l)=\frac{\left\langle A_{n+l} A_n\right\rangle\left\langle A_n\right\rangle^2}{\left\langle A_n{ }^2\right\rangle-\left\langle A_n\right\rangle^2} \\
# \text { where } \left\langle A_{n+l} A_n\right\rangle=\frac{1}{M} \sum_{n=1}^M A_{n+1} A_n .
# \end{array}
import numpy as np
import matplotlib.pyplot as plt

# Autocorrelation function calculation
def autocorrelation(sample_points, lag):
    M = len(sample_points)
    if lag >= M:
        return 0.0

    mean_sample = np.mean(sample_points)
    numerator = sum((sample_points[i] * sample_points[i + lag] for i in range(M - lag))) / M
    denominator = np.var(sample_points)
    
    if denominator == 0:
        return 0.0
    else:
        return (numerator - mean_sample * mean_sample) / denominator

# Calculate autocorrelation values for various lags
max_lag = 89000  
autocorrelation_values = [autocorrelation(total_points, lag) for lag in range(max_lag)]

# Plot the autocorrelation function
plt.figure(figsize=(10, 5))
plt.plot(autocorrelation_values)
plt.title("Autocorrelation Function")
plt.xlabel("Lag")
plt.ylabel("Autocorrelation")
plt.grid(True)
plt.show()



