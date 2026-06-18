#ifndef NUMERICALLIB_NONLINEAR_EQUATIONS_H
#define NUMERICALLIB_NONLINEAR_EQUATIONS_H

#include <functional>
#include <string>
#include <vector>
#include <utility>

namespace NumericalLibrary {

/**
 * @brief Newton's method using an analytical derivative function.
 *
 * @param df_analytical Analytical derivative function f'(x).
 * @param f Function f(x) whose root is sought.
 * @param x0 Initial guess.
 * @param tol Convergence tolerance on successive iterates.
 * @param max_iter Maximum number of iterations.
 * @param outFile Optional output CSV filename for iteration history (can be empty).
 * @return Approximated root.
 *
 * @example
 * double root = newton_analytical(df, f, 1.0, 1e-12, 100, "");
 */
double newton_analytical(std::function<double(double)> df_analytical, std::function<double(double)> f, double x0, double tol, int max_iter, const std::string& outFile);

/**
 * @brief Newton's method where derivative is approximated numerically.
 */
double newton_numerical(std::function<double(double)> f, double x0, double tol, int max_iter, const std::string& outFile);

/**
 * @brief Secant method for root finding requiring two initial guesses.
 */
double secantMethod(std::function<double(double)> f, double x0, double x1, double tolerance, int maxIterations, const std::string& outFile);

/**
 * @brief Numerical derivative using central finite differences.
 *
 * @param f Function to differentiate.
 * @param x Point at which derivative is approximated.
 * @param h Finite-difference step (default unspecified here; see implementation).
 * @return Approximate derivative f'(x).
 */
double derivativeNumerical(const std::function<double(double)>& f, double x, double h);

// Regula falsi (false position) and helpers

/**
 * @brief Regula falsi (false position) method.
 *
 * Construct the secant between (a,f(a)) and (b,f(b)) and use its intersection
 * with the x-axis as the next approximation. This routine requires that f(a)
 * and f(b) have opposite signs (a bracketing interval). The method keeps the
 * bracketing property and returns the final approximation in `root`.
 *
 * @param f Function f(x) whose root is sought.
 * @param a Left end of the bracketing interval.
 * @param b Right end of the bracketing interval.
 * @param root Output parameter set to the computed approximation.
 * @param tol Convergence tolerance (on function value or interval size).
 * @param maxIter Maximum number of iterations.
 * @param verbose If true, prints per-iteration diagnostics to stdout.
 * @param iter_log Optional pointer to a vector to collect (iter,approx) pairs.
 * @return true if an approximation was produced (or converged), false if the
 *         preconditions (sign change) are not met.
 */
bool false_position(const std::function<double(double)>& f, double a, double b, double& root,
	double tol, int maxIter, bool verbose,
	std::vector<std::pair<int,double>>* iter_log = nullptr);

/**
 * @brief Find roots by scanning subintervals using false position.
 *
 * The interval [a,b] is split into l_podpodzialow subintervals. For each
 * subinterval where a sign change (or near-zero sample) is detected, the
 * false_position routine is applied. Found roots are deduplicated and
 * returned in `roots`. Optionally per-root iteration logs are returned via
 * all_iter_logs.
 *
 * @param f Function f(x) whose roots are sought.
 * @param a Left end of the search interval.
 * @param b Right end of the search interval.
 * @param l_podpodzialow Number of subdivisions to scan.
 * @param roots Output vector of distinct root approximations.
 * @param verbose If true, print diagnostic information during scanning.
 * @param tol Convergence / sampling tolerance.
 * @param maxIter Maximum iterations for false_position on each subinterval.
 * @param all_iter_logs Optional pointer to collect per-root iteration logs.
 * @return true if at least one root was found, false otherwise.
 */
bool find_roots_false_position(const std::function<double(double)>& f, double a, double b, int l_podpodzialow,
	std::vector<double>& roots, bool verbose, double tol, int maxIter,
	std::vector<std::vector<std::pair<int,double>>>* all_iter_logs = nullptr);

// Bisection and helper that finds roots by subdividing an interval
/**
 * @brief Bisection (interval halving) root-finder.
 *
 * Given a bracketing interval [a,b] with f(a)*f(b) <= 0, the routine
 * repeatedly halves the interval and keeps the subinterval that contains a
 * sign change. Bisection is robust and guarantees convergence (linear rate).
 *
 * @param f Function f(x) whose root is sought.
 * @param a Left end of the bracketing interval.
 * @param b Right end of the bracketing interval.
 * @param root Output parameter set to the computed approximation.
 * @param tol Convergence tolerance on interval size or function value.
 * @param maxIter Maximum number of iterations.
 * @param verbose If true, prints per-iteration diagnostics.
 * @param iter_log Optional pointer to collect (iter,approx) pairs.
 * @return true if an approximation was produced (or converged), false if the
 *         preconditions (sign change) are not met.
 */
bool bisection(const std::function<double(double)>& f, double a, double b, double& root,
	double tol, int maxIter, bool verbose,
	std::vector<std::pair<int,double>>* iter_log = nullptr);

/**
 * @brief Find roots by subdividing [a,b] and applying bisection on sign-changing subintervals.
 *
 * The function splits the interval into l_podpodzialow parts, checks for
 * sign changes (and near-zero samples) and runs bisection where appropriate.
 * Distinct roots are deduplicated and returned in `roots`. Per-root logs can
 * be returned via all_iter_logs.
 *
 * @param f Function f(x) whose roots are sought.
 * @param a Left end of the search interval.
 * @param b Right end of the search interval.
 * @param l_podpodzialow Number of subdivisions to scan.
 * @param roots Output vector of distinct root approximations.
 * @param verbose If true, prints diagnostics during scanning.
 * @param tol Convergence / sampling tolerance.
 * @param maxIter Maximum iterations for bisection on each subinterval.
 * @param all_iter_logs Optional pointer to collect per-root iteration logs.
 * @return true if at least one root was found, false otherwise.
 */
bool find_roots_bisection(const std::function<double(double)>& f, double a, double b, int l_podpodzialow,
	std::vector<double>& roots, bool verbose, double tol, int maxIter,
	std::vector<std::vector<std::pair<int,double>>>* all_iter_logs = nullptr);

} // namespace NumericalLibrary

#endif // NUMERICALLIB_NONLINEAR_EQUATIONS_H
