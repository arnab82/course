#include <iostream>
#include <vector>
#include <cmath>

// Function for nth order polynomial interpolation
double nth_order_interpolation(double x, std::vector<double> x_data, std::vector<double> y_data) {
    int n = x_data.size() - 1;
    double interpolated_value = 0.0;
    
    for (int i = 0; i <= n; ++i) {
        double term = y_data[i];
        
        for (int j = 0; j <= n; ++j) {
            if (i != j) {
                term *= (x - x_data[j]) / (x_data[i] - x_data[j]);
            }
        }
        
        interpolated_value += term;
    }
    
    return interpolated_value;
}

// Function to calculate the error for nth order polynomial interpolation
double error_for_nth_polynomial(double x, std::vector<double> x_data, std::vector<double> y_data) {
    int n = x_data.size() - 1;
    double error = 0.0;
    double term = nth_order_interpolation(x, x_data, y_data) / std::tgamma(n + 2);
    
    for (int i = 0; i <= n; ++i) {
        term *= (x - x_data[i]);
    }
    
    error += term;
    return error;
}



// Function for Neville's Lagrange interpolation
std::pair<double, std::vector<std::vector<double>>> neville_lagrange_interpolation(double x, std::vector<double> x_data, std::vector<double> y_data) {
    int n = x_data.size() - 1;
    std::vector<std::vector<double>> function_value(n + 1, std::vector<double>(n + 1, 0.0));

    // Initialize the first column of the function_value matrix with y_data
    for (int i = 0; i <= n; ++i) {
        function_value[i][0] = y_data[i];
    }

    // Neville's interpolation
    for (int i = 1; i <= n; ++i) {
        for (int j = i; j <= n; ++j) {
            double numerator1 = (x - x_data[j - i]) * function_value[j][i - 1];
            double numerator2 = (x - x_data[j]) * function_value[j - 1][i - 1];
            double denominator = x_data[j] - x_data[j - i];
            function_value[j][i] = (numerator1 - numerator2) / denominator;
        }
    }

    double approximated_value = function_value[n][n];
    
    return std::make_pair(approximated_value, function_value);
}

double error_neville_lagrange_interpolation(double x, std::vector<double> x_data, std::vector<double> y_data) {
    auto result = neville_lagrange_interpolation(x, x_data, y_data);
    int n = x_data.size() - 1;

    double numerical_uncertainty = 0.5 * (std::abs(result.second[n][n] - result.second[n - 1][n - 1]) +
                                          std::abs(result.second[n][n] - result.second[n][n - 1]));

    return numerical_uncertainty;
}
// Function for linear interpolation
double linear_interpolation(std::vector<double> x_values, std::vector<double> y_values, double x) {
    // Number of data points
    int n = x_values.size();
    
    // Finding the interval which contains x
    int i = 0;
    while (i < n) {
        if (x < x_values[i]) {
            break;
        }
        i++;
    }
    
    // If x is outside the data range, return 0.0
    if (i == 0 || i == n) {
        std::cout << "x is outside the data range and hence the function will show 0 as interpolated value" << std::endl;
        return 0.0;
    }
    
    // Linear interpolation
    double y = y_values[i - 1] + (y_values[i] - y_values[i - 1]) * (x - x_values[i - 1]) / (x_values[i] - x_values[i - 1]);
    return y;
}

int main() {
    std::vector<double> x_datas = {0.2, 0.3, 0.4, 0.5, 0.6};
    std::vector<double> y_datas = {0.0015681, 0.0077382, 0.023579, 0.054849, 0.10696};
    double x1 = 0.16;
    double x2 = 0.46;
    double result1 = nth_order_interpolation(x1, x_datas, y_datas);
    double error1 = error_for_nth_polynomial(x1, x_datas, y_datas);
    double result2 = nth_order_interpolation(x2, x_datas, y_datas);
    double error2 = error_for_nth_polynomial(x2, x_datas, y_datas);

    std::cout << "Polynomial value at x = " << x1 << ": " << result1 << " Error: " << error1 << std::endl;
    std::cout << "Polynomial value at x = " << x2 << ": " << result2 << " Error: " << error2 << std::endl;

    auto result1_ = neville_lagrange_interpolation(x1, x_datas, y_datas);
    std::cout << "Approximated value at x = " << x1 << ": " << result1_.first << std::endl;
    std::cout << "Neville iterations table:" << std::endl;
    for (int i = 0; i < result1_.second.size(); ++i) {
        for (int j = 0; j < result1_.second[i].size(); ++j) {
            std::cout << result1_.second[i][j] << "\t";
        }
        std::cout << std::endl;
    }

    auto result2_ = neville_lagrange_interpolation(x2, x_datas, y_datas);
    std::cout << "Approximated value at x = " << x2 << ": " << result2_.first << std::endl;
    std::cout << "Neville iterations table:" << std::endl;
    for (int i = 0; i < result2_.second.size(); ++i) {
        for (int j = 0; j < result2_.second[i].size(); ++j) {
            std::cout << result2_.second[i][j] << "\t";
        }
        std::cout << std::endl;
    }

    double uncertainty1 = error_neville_lagrange_interpolation(x1, x_datas, y_datas);
    double uncertainty2 = error_neville_lagrange_interpolation(x2, x_datas, y_datas);
    std::cout << "Numerical uncertainty at x = " << x1 << ": " << uncertainty1 << std::endl;
    std::cout << "Numerical uncertainty at x = " << x2 << ": " << uncertainty2 << std::endl;
    double interpolated_value1 = linear_interpolation(x_datas, y_datas, x1);
    double interpolated_value2 = linear_interpolation(x_datas, y_datas, x2);

    std::cout << "Interpolated value at x = " << x1 << ": " << interpolated_value1 << std::endl;
    std::cout << "Interpolated value at x = " << x2 << ": " << interpolated_value2 << std::endl;


    return 0;
}
