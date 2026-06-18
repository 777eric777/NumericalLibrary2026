#ifndef NUMERICALLIB_APPROXIMATION_H
#define NUMERICALLIB_APPROXIMATION_H

#include <vector>
#include <functional>

namespace NumericalLibrary {

/**
 * @brief Compute least-squares polynomial approximation coefficients for a
 * given function f on interval [range_a, range_b]. Uses Gaussian quadrature
 * to compute integrals required for normal equations.
 *
 * @param f Function to approximate.
 * @param range_a Left endpoint of approximation interval.
 * @param range_b Right endpoint of approximation interval.
 * @param degree Degree of polynomial approximation.
 * @param gl_n Number of Gauss-Legendre nodes for inner quadrature (default 8).
 * @param gl_k Number of subintervals for composite Gauss-Legendre (default 10).
 * @return Vector of coefficients [a0, a1, ..., adegree].
 *
 * @example
 * auto coeffs = leastSquaresApprox([](double x){ return sin(x); }, 0.0, 3.14, 5);
 */
std::vector<double> leastSquaresApprox(const std::function<double(double)>& f, double range_a, double range_b, int degree, int gl_n = 8, int gl_k = 10);

} // namespace NumericalLibrary

#endif // NUMERICALLIB_APPROXIMATION_H
