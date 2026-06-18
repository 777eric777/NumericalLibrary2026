#ifndef NUMERICALLIB_INTERPOLATION_H
#define NUMERICALLIB_INTERPOLATION_H

#include <vector>

namespace NumericalLibrary {

/**
 * @brief Evaluate the Lagrange interpolation polynomial at x given nodes xi
 * and function values fxi.
 *
 * @param x Point where the polynomial is evaluated.
 * @param xi Vector of interpolation nodes.
 * @param fxi Vector of function values at the nodes.
 * @return Interpolated value P(x).
 *
 * @example
 * double y = lagrange(1.23, xi, fxi);
 */
double lagrange(double x, const std::vector<double>& xi, const std::vector<double>& fxi);

/**
 * @brief Compute Newton divided-differences coefficients for nodes xi and
 * values fxi. Returns vector of coefficients for Newton form.
 *
 * @param xi Nodes.
 * @param fxi Function values at nodes.
 * @param size Number of nodes.
 * @return Vector of divided-difference coefficients.
 *
 * @example
 * auto coeffs = i_newtona(xi, fxi, xi.size());
 */
std::vector<double> i_newtona(const std::vector<double>& xi, const std::vector<double>& fxi, int size);

/**
 * @brief Evaluate polynomial using Horner's method (coefficients in
 * increasing order).
 *
 * Implemented as an inline function to avoid duplicate symbol definitions when
 * the header is included from multiple translation units.
 */
inline double p_hornera(double x, const std::vector<double>& wspolczynniki) {
	double wynik = 0.0;
	for (int i = (int)wspolczynniki.size() - 1; i >= 0; --i) wynik = wynik * x + wspolczynniki[i];
	return wynik;
}

/**
 * @brief Evaluate polynomial in the natural basis (sum ai * x^i).
 */
inline double p_naturalna(double x, const std::vector<double>& wspolczynniki) {
	double wynik = 0.0;
	for (int i = 0; i < (int)wspolczynniki.size(); ++i) wynik += wspolczynniki[i] * pow(x, i);
	return wynik;
}

} // namespace NumericalLibrary

#endif // NUMERICALLIB_INTERPOLATION_H
