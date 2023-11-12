import random
import math
import numpy as np

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


def monte_carlo_integration(a, b, M,N,skipped_samples):
    #initialize the variables
    integral_sum = 0.0
    sample_sum = 0.0
    total_weight = 0.0
    autocorrelation_values = []
    sample_squared_sum = 0.0  
    count=0
    count1=0
    count2=0

    for i in range(M + skipped_samples):
        
        count+=1
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
            count1+=1
        # print(x)
        integrand_val = f(x)
        autocorrelation_values.append(integrand_val/(W_normalized(x,N)))
        if i >= 80000:
            integral_sum += integrand_val /(W_normalized(x,N))
            sample_sum += integrand_val/(W_normalized(x,N))
            sample_squared_sum += integrand_val**2
            total_weight += (W_normalized(x,N))
            count2+=1
            
    acceptance_ratio=count2/count
    
    
    integral_estimate = integral_sum /M

    # Calculate the sample mean and the standard deviation of the sample mean
    sample_mean = sample_sum / M
    sample_mean_squared = sample_squared_sum / M
    sample_mean_stddev = math.sqrt((sample_mean_squared - sample_mean**2) / M)

    return integral_estimate, sample_mean, sample_mean_stddev,autocorrelation_values,acceptance_ratio

a = 0.0  # Lower limit of the integration interval
b = 1.0  # Upper limit of the integration interval
M = 100000  # Number of sample points
N = estimate_normalization_constant(M, a, b)

print("Normalization Constant:", N)
estimated_integral, mean_sample, sample_mean_stddev,autocorrelation_values ,acceptance_ratio= monte_carlo_integration(a, b, M,N,80000)
print("the accepatnce rate of importannt sampling is ",acceptance_ratio*100,"%")
print("Monte Carlo Integral Estimate with Importance Sampling:", estimated_integral)
print("Sample Mean:", mean_sample)
print("Standard Deviation of Sample Mean:", sample_mean_stddev)
for i in (50000,55000,60000,65000,70000,75000,80000,85000,90000,95000,100000):
    values=monte_carlo_integration(a, b, M,N,i)
    print("the integral is for skipped samples ",i,"is ",values[0])
print("the integral values are changing with the skipped samples\n",
      "because the delta value for selecting x_next should also change with the skipped samples\n"
      "hence we should adjust the delta value according to the skipped samples\n ")
