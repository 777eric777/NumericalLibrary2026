#ifndef NUMERICALLIB_ODE_H
#define NUMERICALLIB_ODE_H

#include <vector>
#include <functional>

namespace NumericalLibrary {

/**
 * @brief Explicit Euler integration for the scalar ODE y' = f(t, y).
 *
 * @param f Right-hand side function f(t, y).
 * @param t0 Initial time.
 * @param y0 Initial state y(t0).
 * @param h Step size.
 * @param t_end End time.
 * @return Pair of vectors (time, values) containing discrete solution points.
 *
 * @example
 * auto [t, y] = eulerExplicit(f, 0.0, 1.0, 0.01, 1.0);
 */
std::pair<std::vector<double>, std::vector<double>> eulerExplicit(const std::function<double(double, double)>& f, double t0, double y0, double h, double t_end);

/**
 * @brief Heun's (improved Euler) method for scalar ODEs.
 */
std::pair<std::vector<double>, std::vector<double>> heunExplicit(const std::function<double(double, double)>& f, double t0, double y0, double h, double t_end);

/**
 * @brief Midpoint method (second-order Runge-Kutta).
 */
std::pair<std::vector<double>, std::vector<double>> midpointMethod(const std::function<double(double, double)>& f, double t0, double y0, double h, double t_end);

/**
 * @brief Classical 4th-order Runge-Kutta integrator for scalar ODEs.
 */
std::pair<std::vector<double>, std::vector<double>> rungeKutta4(double y0, double t0, double t_end, double h, const std::function<double(double,double)>& f);

} // namespace NumericalLibrary

#endif // NUMERICALLIB_ODE_H
