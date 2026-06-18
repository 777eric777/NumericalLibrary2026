#ifndef NUMERICALLIB_LINEAR_SYSTEMS_H
#define NUMERICALLIB_LINEAR_SYSTEMS_H

#include <vector>

namespace NumericalLibrary {

/**
 * @brief Solve linear systems and related utilities
 */

/**
 * @brief Performs Gaussian elimination (with partial pivoting) to transform
 * the linear system Ax = b into an upper-triangular system.
 *
 * @param N The dimension of the system (number of equations).
 * @param a The coefficient matrix (will be modified by the routine).
 * @param b The right-hand side vector (will be modified accordingly).
 * @return The augmented matrix after elimination (for debugging/inspection).
 *
 * @example
 * std::vector<std::vector<double>> A = { ... };
 * std::vector<double> b = { ... };
 * auto result = gauss(N, A, b);
 */
std::vector<std::vector<double>> gauss(int N, std::vector<std::vector<double>>& a, std::vector<double>& b);

/**
 * @brief Back-substitution to recover the solution vector from an
 * upper-triangular system produced by Gaussian elimination.
 *
 * @param a Upper-triangular matrix.
 * @param b Right-hand side vector.
 * @param n Dimension of the system.
 * @return Solution vector x such that a*x = b.
 *
 * @example
 * auto x = xGauss(upperA, b, n);
 */
std::vector<double> xGauss(const std::vector<std::vector<double>>& a, const std::vector<double>& b, int n);

/**
 * @brief Compute LU decomposition of matrix A such that A = L * U.
 *
 * @param A Input square matrix of size N x N.
 * @param N Dimension of the matrix.
 * @return Pair of matrices (L, U) where L is lower-triangular with unit diagonal
 * and U is upper-triangular.
 *
 * @example
 * auto [L, U] = LU_decomposition(A, N);
 */
std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> LU_decomposition(const std::vector<std::vector<double>>& A, int N);

/**
 * @brief Forward substitution to solve L*y = b for y when L is lower-triangular.
 *
 * @param L Lower-triangular matrix with unit or non-unit diagonal.
 * @param b Right-hand side vector.
 * @param n Dimension of the system.
 * @return Vector y solution of L*y = b.
 *
 * @example
 * auto y = forwardSub(L, b, n);
 */
std::vector<double> forwardSub(const std::vector<std::vector<double>>& L, const std::vector<double>& b, int n);

/**
 * @brief Compute the residual norm ||Ax - b||_2 for a computed solution x.
 *
 * @param A Coefficient matrix.
 * @param x Solution vector.
 * @param b Right-hand side vector.
 * @return Euclidean norm of the residual vector.
 *
 * @example
 * double res = computeResidualNorm(A, x, b);
 */
double computeResidualNorm(const std::vector<std::vector<double>>& A, const std::vector<double>& x, const std::vector<double>& b);

} // namespace NumericalLibrary

#endif // NUMERICALLIB_LINEAR_SYSTEMS_H
