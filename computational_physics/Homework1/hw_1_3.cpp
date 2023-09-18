#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

// Function to calculate f(x)
double f(double x) {
    return exp(-x);
}

// Analytical integral
double analytical_integral(double a, double b) {
    return exp(-a) - exp(-b);
}

// Simpson's rule for numerical integration
double simpsons_rule(double (*f)(double), double a, double b, int n) {
    double h = (b - a) / n;
    std::vector<double> x_values(n + 1);

    for (int i = 0; i <= n; ++i) {
        x_values[i] = a + i * h;
    }

    double integral = 0.0;
    
    for (int i = 1; i <= n / 2; ++i) {
        integral += (h / 3) * (f(x_values[2 * i - 1]) + 4 * f(x_values[2 * i]) + f(x_values[2 * i + 1]));
    }

    if (n % 2 == 0) {
        integral += 0.0;
    } else {
        integral += (h / 12) * (-f(x_values[n - 1]) + 8 * f(x_values[n]) + 5 * f(x_values[n + 1]));
    }

    return integral;
}

int main() {
    double a = 0.0;
    double b = 1.0;
    double exact_integral = analytical_integral(a, b);
    
    int n_list[] = {998, 999, 1000, 1001, 10000, 9999};
    int random_index = rand() % (sizeof(n_list) / sizeof(n_list[0]));
    int n = n_list[random_index];
    std::cout << "Randomly generated n: " << n << std::endl;
    
    double numerical_integral = simpsons_rule(f, a, b, n);
    double numerical_uncertainty = std::abs(exact_integral - numerical_integral);

    std::cout << "Analytical Integral: " << exact_integral << std::endl;
    std::cout << "Numerical Integral (Simpson's Rule): " << numerical_integral << std::endl;
    std::cout << "Numerical uncertainty: " << numerical_uncertainty << std::endl;

    printf("n = %4.2d Analytical Integral = %16.12f numerical integral = %16.12f numerical uncertainty = %16.15f\n",
           n, exact_integral, numerical_integral, numerical_uncertainty);

    return 0;
}
