#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

// Function for partial pivoting
void partial_pivoting(std::vector<std::vector<double>>& A, std::vector<double>& b, int col) {
    int n = A.size();
    int max_row = col;
    double max_val = std::abs(A[col][col]);

    // Finding the row with the maximum value in the current column
    for (int i = col + 1; i < n; ++i) {
        double val = std::abs(A[i][col]);
        if (val > max_val) {
            max_row = i;
            max_val = val;
        }
    }

    // Swap rows in A and b
    std::swap(A[col], A[max_row]);
    std::swap(b[col], b[max_row]);
}

// Function for Gaussian elimination
std::vector<double> gauss_elimination(std::vector<std::vector<double>> A, std::vector<double> b) {
    int n = A.size();

    // Forward elimination
    for (int i = 0; i < n; ++i) {
        partial_pivoting(A, b, i);
        for (int j = i + 1; j < n; ++j) {
            double factor = A[j][i] / A[i][i];
            for (int k = i; k < n; ++k) {
                A[j][k] -= factor * A[i][k];
            }
            b[j] -= factor * b[i];
        }
    }

    // Backward substitution
    std::vector<double> x(n);
    for (int i = n - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < n; ++j) {
            sum += A[i][j] * x[j];
        }
        x[i] = (b[i] - sum) / A[i][i];
    }

    return x;
}

int main() {
    // Given matrix A
    std::vector<std::vector<double>> A = {
        {1.0, -1.0, -2.0, 1.0},
        {3.0, 2.0, -1.0, 2.0},
        {2.0, 3.0, 1.0, -3.0},
        {10.0, -4.0, 3.0, 2.0}
    };

    // Right-hand side vector b
    std::vector<double> b = {1.0, 4.0, 2.0, 3.0};

    // Solving the system of equations using Gaussian elimination with partial pivoting
    std::vector<double> x = gauss_elimination(A, b);

    // Printing the solution
    std::cout << "Evaluated Solution (x, y, z, u): ";
    for (size_t i = 0; i < x.size(); ++i) {
        std::cout << std::fixed << std::setprecision(4) << x[i] << " ";
    }
    std::cout << std::endl;

    // Expected solution
    std::vector<double> expected_solution = {21.0 / 34.0, 43.0 / 68.0, -13.0 / 34.0, 1.0 / 4.0};
    std::cout << "Expected Solution (x, y, z, u): ";
    for (size_t i = 0; i < expected_solution.size(); ++i) {
        std::cout << std::fixed << std::setprecision(4) << expected_solution[i] << " ";
    }
    std::cout << std::endl;

    return 0;
}
