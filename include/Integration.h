#ifndef NUMERICALLIB_INTEGRATION_H
#define NUMERICALLIB_INTEGRATION_H

#include <vector>
#include <functional>

namespace NumericalLibrary {

/**
 * @brief Composite trapezoidal integration for a polynomial represented by
 * coefficients ai over the interval defined by range (size 2: [a,b]).
 *
 * @param ai Polynomial coefficients (a0, a1, ...).
 * @param range Two-element vector with [a, b].
 * @param liczbaPodprzedzialow Number of subintervals (default 1000).
 * @return Approximated integral value.
 */
double integralTrapez(const std::vector<double>& ai, const std::vector<double>& range, int liczbaPodprzedzialow = 1000);

/**
 * @brief Composite Simpson integration for a polynomial represented by
 * coefficients ai over interval [a,b].
 */
double integralSimpsona(const std::vector<double>& ai, const std::vector<double>& range, int liczbaPodprzedzialow = 1000);

/**
 * @brief Composite trapezoidal integration for an arbitrary function f on [a,b].
 *
 * @param f Function to integrate.
 * @param a Left endpoint.
 * @param b Right endpoint.
 * @param liczbaPodprzedzialow Number of subintervals (default 1000).
 * @return Approximated integral value.
 */
double integralTrapez(const std::function<double(double)>& f, double a, double b, int liczbaPodprzedzialow = 1000);

/**
 * @brief Composite Simpson integration for an arbitrary function f on [a,b].
 */
double integralSimpsona(const std::function<double(double)>& f, double a, double b, int liczbaPodprzedzialow = 1000);

/**
 * @brief Gauss-Legendre quadrature with n nodes on interval [a,b].
 *
 * @param f Function to integrate.
 * @param a Left endpoint.
 * @param b Right endpoint.
 * @param n Number of Gauss nodes.
 * @return Approximated integral value.
 */
double gaussLegendre(const std::function<double(double)>& f, double a, double b, int n);

/**
 * @brief Composite Gauss-Legendre over k subintervals using n nodes per subinterval.
 */
double gaussLegendreComposite(const std::function<double(double)>& f, double a, double b, int n, int k);

} // namespace NumericalLibrary

#endif // NUMERICALLIB_INTEGRATION_H
