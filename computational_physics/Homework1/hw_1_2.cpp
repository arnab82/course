#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

// Function to calculate f(x)
double f(double x) {
    return x * x * x * cos(x);
}

// First-order derivative using three-point formula
double first_derivative(double x, double h) {
    return (f(x + h) - f(x - h)) / (2 * h);
}

// Second-order derivative using three-point formula
double second_derivative(double x, double h) {
    return (f(x + h) - 2 * f(x) + f(x - h)) / (h * h);
}

// Analytical first derivative
double analytical_first_derivative(double x) {
    return 3 * x * x * cos(x) - x * x * x * sin(x);
}

// Analytical second derivative
double analytical_second_derivative(double x) {
    return 6 * x * cos(x) - 6 * x * x * sin(x) - x * x * x * cos(x);
}

int main() {
    double x_start = M_PI / 2;
    double x_end = M_PI;
    int num_intervals = 100;
    double h = (x_end - x_start) / num_intervals;

    std::vector<double> x_values;
    std::vector<double> numerical_first_derivatives;
    std::vector<double> numerical_second_derivatives;
    std::vector<double> analytical_first_derivatives;
    std::vector<double> analytical_second_derivatives;
    std::vector<double> numerical_errorfirst;
    std::vector<double> numerical_errorsecond;

    for (int i = 0; i <= num_intervals; ++i) {
        double x = x_start + i * h;

        // Calculating derivatives using three-point formulas
        double first_deriv = first_derivative(x, h);
        double second_deriv = second_derivative(x, h);

        // Calculating analytical derivatives
        double analytical_first = analytical_first_derivative(x);
        double analytical_second = analytical_second_derivative(x);

        // Calculating numerical accuracy
        double numerical_error_first = std::abs(analytical_first - first_deriv);
        double numerical_error_second = std::abs(analytical_second - second_deriv);

        x_values.push_back(x);
        numerical_first_derivatives.push_back(first_deriv);
        numerical_second_derivatives.push_back(second_deriv);
        analytical_first_derivatives.push_back(analytical_first);
        analytical_second_derivatives.push_back(analytical_second);
        numerical_errorfirst.push_back(numerical_error_first);
        numerical_errorsecond.push_back(numerical_error_second);
    }

    std::cout << "x\t\tNumerical 1st Derivative\tNumerical 2nd Derivative\tNumerical Error (1st)\tNumerical Error (2nd)" << std::endl;
    for (size_t i = 0; i < x_values.size(); ++i) {
        std::cout << std::fixed << std::setprecision(6)
                  << x_values[i] << "\t\t" << numerical_first_derivatives[i] << "\t\t" << numerical_second_derivatives[i]
                  << "\t\t\t\t" << numerical_errorfirst[i] << "\t\t" << numerical_errorsecond[i] << std::endl;
    }

    return 0;
}
