using Random
using Plots
using Statistics
function f(x)
    return x^2
end

function W(x)
    return x^2 + 2.0
end

function W_normalized(x, N)
    return W(x) / N
end

function estimate_normalization_constant(M, a, b)
    total_weight = 0.0
    for _ in 1:M
        x = a + (b - a) * rand() 
        weight = W(x)
        total_weight += weight
    end
    N = (b - a) * total_weight / M
    return N
end

function monte_carlo_integration(a, b, M, N)
    integral_sum = 0.0
    sample_sum = 0.0
    total_weight = 0.0
    autocorrelation_values = Float64[]
    sample_squared_sum = 0.0

    for i in 1:(M + 80000)
        x = rand()  # Generate a random number in [0, 1]
        # print(x)
        weight_0 = W_normalized(x, N)
        rand_val = (2*rand()-1.0) *1.27
        x_next = (b - a) * rand_val + x
        weight_1 = W_normalized(x_next, N)
        prob = weight_1 / weight_0
        si = rand()
        
        if prob >= si && x_next <= 1.0
            x = x_next
        end

        integrand_val = f(x)
        push!(autocorrelation_values, integrand_val / W_normalized(x, N))

        if i >= 80000
            integral_sum += integrand_val / W_normalized(x, N)
            sample_sum += integrand_val / W_normalized(x, N)
            sample_squared_sum += integrand_val^2
            total_weight += W_normalized(x, N)
        end
    end

    integral_estimate = integral_sum / M
    sample_mean = sample_sum / M
    sample_mean_squared = sample_squared_sum / M
    sample_mean_stddev = sqrt((sample_mean_squared - sample_mean^2) / M)

    return integral_estimate, sample_mean, sample_mean_stddev, autocorrelation_values
end

function autocorrelation(sample_points, lag)
    M = length(sample_points)
    if lag >= M
        return 0.0
    end

    mean_sample = mean(sample_points)
    numerator = sum([sample_points[i] * sample_points[i + lag] for i in 1:(M - lag)]) / M
    denominator = var(sample_points)

    if denominator == 0.0
        return 0.0
    else
        return (numerator - mean_sample^2) / denominator
    end
end

a = 0.0  # Lower limit of the integration interval
b = 1.0  # Upper limit of the integration interval
M = 100000  # Number of sample points
N = estimate_normalization_constant(M, a, b)
println("Normalization Constant: $N")
estimated_integral, mean_sample, sample_mean_stddev, autocorrelation_values = monte_carlo_integration(a, b, M, N)
println("Monte Carlo Integral Estimate with Importance Sampling: $estimated_integral")
println("Sample Mean: $mean_sample")
println("Standard Deviation of Sample Mean: $sample_mean_stddev")

max_lag = 89000  
autocorrelation_l = [autocorrelation(autocorrelation_values, lag) for lag in 1:max_lag]

# Plot the autocorrelation function
plot(autocorrelation_l, xlabel="Lag", ylabel="Autocorrelation", title="Autocorrelation Function", legend=false, grid=true)
