#include <iostream>
#include <vector>

// Function to perform LU decomposition with pivoting
void lu_decomposition_pivoting(std::vector<std::vector<double>>& A, std::vector<int>& P) {
    int n = A.size();
    P.resize(n);

    for (int i = 0; i < n; i++) {
        P[i] = i;
    }

    for (int k = 0; k < n - 1; k++) {
        int pivot_row = k;

        for (int i = k + 1; i < n; i++) {
            if (std::abs(A[i][k]) > std::abs(A[pivot_row][k])) {
                pivot_row = i;
            }
        }

        if (pivot_row != k) {
            std::swap(P[k], P[pivot_row]);
            std::swap(A[k], A[pivot_row]);
        }

        for (int i = k + 1; i < n; i++) {
            double factor = A[i][k] / A[k][k];
            A[i][k] = factor;

            for (int j = k + 1; j < n; j++) {
                A[i][j] -= factor * A[k][j];
            }
        }
    }
}

// Function to solve a linear system using LU decomposition with pivoting
std::vector<double> solve_linear_system(const std::vector<std::vector<double>>& A, const std::vector<int>& P, const std::vector<double>& b) {
    int n = A.size();
    std::vector<double> x(n, 0.0);
    std::vector<double> y(n);

    // Forward substitution (Ly = Pb)
    for (int i = 0; i < n; i++) {
        double sum = 0.0;

        for (int j = 0; j < i; j++) {
            sum += A[i][j] * y[j];
        }

        y[i] = b[P[i]] - sum;
    }

    // Backward substitution (Ux = y)
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;

        for (int j = i + 1; j < n; j++) {
            sum += A[i][j] * x[j];
        }

        x[i] = (y[i] - sum) / A[i][i];
    }

    return x;
}

int main() {
    // Given matrix A and right-hand side vector b
    std::vector<std::vector<double>> A = {{1.0, -1.0, -2.0, 1.0},
                                          {3.0, 2.0, -1.0, 2.0},
                                          {2.0, 3.0, 1.0, -3.0},
                                          {10.0, -4.0, 3.0, 2.0}};

    std::vector<double> b = {1.0, 4.0, 2.0, 3.0};
    std::vector<int> permutation;
    
    // Perform LU decomposition with pivoting
    lu_decomposition_pivoting(A, permutation);

    // Solve the linear system
    std::vector<double> x = solve_linear_system(A, permutation, b);

    // Display the solution vector
    std::cout << "Solution vector:" << std::endl;
    for (int i = 0; i < x.size(); i++) {
        std::cout << "x[" << i << "] = " << x[i] << std::endl;
    }

    return 0;
}
