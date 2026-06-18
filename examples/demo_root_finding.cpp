#include "NumericalLibrary.h"
#include <iostream>

// Demo: root-finding utilities
// This example shows practical usage of three root-finding methods provided by
// the library: newton_analytical, newton_numerical and secantMethod.
//
// Function: f(x) = x^2 - 2 (root at sqrt(2))
// Arguments explained:
// - df (for newton_analytical): analytical derivative f'(x) = 2x
// - x0 / x1: initial guesses for iterative methods
// - tol: convergence tolerance for stopping criterion
// - max_iter: maximum allowed iterations
// - outFile: CSV filename for iteration history (pass empty string to skip file output)

int main() {
	using namespace NumericalLibrary;

	auto f = [](double x) { return x * x - 2.0; };
	auto df = [](double x) { return 2.0 * x; };

	double tol = 1e-12;
	int max_iter = 50;

	// Newton using analytical derivative:
	// Good when you can provide exact derivative; usually converges quadratically
	double root1 = newton_analytical(df, f, 1.0, tol, max_iter, "");
	std::cout << "newton_analytical root = " << root1 << "\n";

	// Newton using numerical derivative:
	// Useful when derivative is hard to compute; less accurate and may be slower
	double root2 = newton_numerical(f, 1.0, tol, max_iter, "");
	std::cout << "newton_numerical root = " << root2 << "\n";

	// Secant method:
	// Requires two initial guesses; does not require derivative
	double root3 = secantMethod(f, 1.0, 2.0, tol, max_iter, "");
	std::cout << "secantMethod root = " << root3 << "\n";

	return 0;
}
