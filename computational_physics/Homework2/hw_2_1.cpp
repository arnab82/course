#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include <ctime>

// Function to find roots using the secant method
double secant_method(double (*equation)(double, double), double a, std::vector<double> x_values, double tol, int max_iter) {
    // Randomly select initial points
    int random_index1 = rand() % x_values.size();
    int random_index2 = rand() % x_values.size();
    double x_i = x_values[random_index1];
    double x_j = x_values[random_index2];
    
    if (equation(x_i, a) * equation(x_j, a) < 0.0) {
        for (int i = 1; i <= max_iter; ++i) {
            double f_i = equation(x_i, a);
            double f_j = equation(x_j, a);
            
            if (std::abs(x_j - x_i) > tol) {
                double x_next = x_j - f_j * (x_j - x_i) / (f_j - f_i);
                x_i = x_j;
                x_j = x_next;
                std::cout << "Iteration No: " << i << " Abs(xj-xi): " << std::fixed << std::setprecision(12) << std::abs(x_j - x_i) << " Xj: " << x_j << std::endl;
            } else if (i == max_iter) {
                std::cout << "It didn't converge" << std::endl;
            } else {
                return x_j;  // Converged to a solution within tolerance
            }
        }
    } else {
        std::cout << "The initial points are not suitable for this method" << std::endl;
    }

    return NAN;  // Return NaN to indicate no root found
}

int main() {
    srand(static_cast<unsigned>(time(0)));
    
    double a = 5.0;
    double tolerance = 1e-8;
    int max_iterations = 100;
    std::vector<double> x_values;

    // Create a range of x values
    double step = 0.1;
    for (double x = M_PI / 2 + tolerance; x <= 3 * M_PI / 2 - tolerance; x += step) {
        x_values.push_back(x);
    }

    // Collect the results in a vector
    std::vector<double> random_results;
    
    for (int i = 1; i <= 20; ++i) {
        double result = secant_method([](double x, double a) { return tan(x) - a / x; }, a, x_values, tolerance, max_iterations);
        random_results.push_back(result);
    }

    std::cout << "\nThe list of approximated roots for 20 random initial guesses:" << std::endl;
    for (size_t i = 0; i < random_results.size(); ++i) {
        if (!std::isnan(random_results[i])) {
            std::cout << "Root " << i + 1 << ": " << random_results[i] << std::endl;
        } else {
            std::cout << "Root " << i + 1 << ": No root found" << std::endl;
        }
    }

    return 0;
}
