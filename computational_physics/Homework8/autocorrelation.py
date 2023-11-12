import random
import math
import numpy as np
import matplotlib.pyplot as plt
def f(x):
    # Define the function to be integrated
    return x**2
def W(x):
    # Define the weight function for importance sampling without normalization
    return x**2 +2.0


def W_normalized(x,N):
    # Define the weight function for importance sampling with normalization
    return W(x)/N

def estimate_normalization_constant(M, a, b):
    total_weight = 0.0
    for _ in range(M):
        x = random.uniform(a, b) 
        weight = W(x)
        total_weight += weight
    N = (b - a) * total_weight / M
    return N

# Calculate the autocorrelation function
def autocorrelation(M, lag,N):
    sample_points = []

    for i in range(0,M+lag):
        x = random.uniform(0, 1)
        weight_0 = W_normalized(x,N)
        rand1 = random.random() 
        rand=(2*rand1-1)*1.27
        x_next = (b - a) * rand + x  # Map the random number to the interval [a, b]
        weight_1 = W_normalized(x_next,N)
        prob=weight_1/weight_0
        si=random.uniform(0,1)
        # print(prob,si,x,x_next)
        
        if prob >=si and x_next<=1.0 :
            x=x_next
        # print(x)
        integrand_val = f(x)
        sample_points.append(integrand_val/(W_normalized(x,N)))    

    mean_sample = np.mean(sample_points[lag:lag+M])
    mean_Sample_l_m=np.mean(sample_points)
    numerator = mean_Sample_l_m*mean_sample-mean_sample**2
    denominator = np.var(sample_points[lag:lag+M])
    if denominator == 0:
        return 0.0
    else:
        # print(numerator/denominator)
        print("no of skipped points", i)
        return numerator / denominator
        
M=10000
a=0.0
b=1.0
N=estimate_normalization_constant(M,0,1)
# Calculate autocorrelation values for various lags
max_lag = 10000
autocor_values = [autocorrelation(M, lag,N) for lag in range(0,max_lag+M)]

# Plot the autocorrelation function
plt.figure(figsize=(10, 5))
plt.plot(autocor_values)
plt.title("Autocorrelation Function")
plt.xlabel("Lag")
plt.ylabel("Autocorrelation")

plt.grid(True)
plt.savefig("C_l_vs_l.png")
plt.show()