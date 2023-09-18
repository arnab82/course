#include <iostream>
#include <vector>

// Function for Gaussian elimination
std::vector<double> gauss_elimination(std::vector<std::vector<double>> A, std::vector<double> b) {
    int n = b.size();
    std::vector<std::vector<double>> Ab(n, std::vector<double>(n + 1));

    // Initialize the augmented matrix Ab
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            Ab[i][j] = A[i][j];
        }
        Ab[i][n] = b[i];
    }

    // Forward elimination
    for (int k = 0; k < n; ++k) {
        // Finding the pivot row and swapping rows if required (pivoting)
        int pivotrow = k;
        for (int i = k + 1; i < n; ++i) {
            if (std::abs(Ab[i][k]) > std::abs(Ab[pivotrow][k])) {
                pivotrow = i;
            }
        }
        std::swap(Ab[k], Ab[pivotrow]);

        // Eliminate the entries below the pivot
        for (int i = k + 1; i < n; ++i) {
            double factor = Ab[i][k] / Ab[k][k];
            for (int j = k; j <= n; ++j) {
                Ab[i][j] -= factor * Ab[k][j];
            }
        }
    }

    // Backward substitution
    std::vector<double> x(n);
    for (int i = n - 1; i >= 0; --i) {
        x[i] = Ab[i][n];
        for (int j = i + 1; j < n; ++j) {
            x[i] -= Ab[i][j] * x[j];
        }
        x[i] /= Ab[i][i];
    }

    return x;
}

int main() {
    // Coefficient matrix A
    std::vector<std::vector<double>> A = {
        {10.0, -4.0, 3.0, 2.0},
        {1.0, -1.0, -2.0, 1.0},
        {3.0, 2.0, -1.0, 2.0},
        {2.0, 3.0, 1.0, -3.0}
    };

    // Right-hand side vector b
    std::vector<double> b = {3.0, 1.0, 4.0, 2.0};

    // Solving the system of equations using Gaussian elimination
    std::vector<double> x = gauss_elimination(A, b);

    // Printing the solution
    std::cout << "Solution:" << std::endl;
    for (size_t i = 0; i < x.size(); ++i) {
        std::cout << "x[" << i << "] = " << x[i] << std::endl;
    }

    return 0;
}
